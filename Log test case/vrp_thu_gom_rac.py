# -*- coding: utf-8 -*-
"""
================================================================================
 TỐI ƯU TUYẾN XE THU GOM RÁC — VYLT 2026
 Cấu trúc: Multi-Trip MD-OVRPTW + VRPMD + Gate Check-in Batching + Service Time
           (with Hybrid Pre-processing)
 Công cụ : Python + Google OR-Tools (định tuyến) + module batching hậu kỳ (C8/C9)
--------------------------------------------------------------------------------
 CÁCH DÙNG:
   pip install ortools matplotlib
   python vrp.py
 File chạy được ngay với DỮ LIỆU MẪU (tổng hợp). Chỗ cắm DỮ LIỆU THẬT được đánh
 dấu bằng "### <<< CẮM DỮ LIỆU THẬT" (ma trận d_ij, t_ij từ Google API; q_i từ ML).
================================================================================
"""
import math, random
from ortools.constraint_solver import pywrapcp, routing_enums_pb2

# ==============================================================================
# 1) THAM SỐ MÔ HÌNH  (ký hiệu theo đúng file đề bài)
# ==============================================================================
class P:
    # --- Khung giờ [20:00 - 24:00], phút tính từ 20:00 ---
    SHIFT_START_MIN   = 0            # 20:00  (a_k >= 20:00)
    SHIFT_END_MIN     = 240          # 24:00
    DT_MAX            = 60           # ΔT_max: cho về trễ tối đa 60' (C2), bị phạt
    HORIZON           = 240 + 60     # 24:00 + ΔT_max

    # --- Xe (xe nâng rác cơ giới) ---
    Q                 = 9000         # Q_k: sức chứa (kg)  ### <<< CẮM DỮ LIỆU THẬT
    THETA             = 0.94         # θ: hệ số an toàn tải (0.93-0.95) -> C1
    K_MAX             = 8            # K_max: fleet tối đa
    SPARE_RATIO       = 0.15         # C10: 10-15% xe standby

    # --- Thời gian phục vụ & cổng KS ---
    S_KS              = 12           # s_KS: cân + đổ tại KS (10-15') (phút)
    DT_KS             = 3            # Δt_KS: dung sai sự cố cổng (phút)
    T_GATE            = 0.5          # t_gate: qua cổng 1 xe = 30s = 0.5'
    B_KS              = 3            # B_KS: tối đa 3 xe / nhóm check-in
    DT_GAP            = 2            # Δt_gap: giãn cách giữa 2 nhóm (phút)
    T_TIMEOUT         = 30          # T_timeout: quá 30' chờ+bãi -> "sự cố" (C9)

    SPEED_KMH         = 25           # tốc độ đêm (đổi d_ij -> t_ij nếu chưa có t thật)

    # --- Trọng số hàm mục tiêu:  min(α·Σd + β·ΣT + γ·K_active + P·Σe_k) ---
    #     (đây là TRỌNG SỐ điều phối, không phải chi phí/km tuyệt đối; tùy chỉnh)
    ALPHA             = 1            # α: theo mét quãng đường
    BETA              = 80           # β: theo phút thời gian
    GAMMA             = 4000         # γ: mỗi xe điều động (chi phí cố định) — tunable
    PENALTY           = 400000       # P: phạt mỗi phút trễ sau 24:00 (khổng lồ)

    RELOAD_COPIES     = 18           # số "bản sao KS" cho phép quay lại đổ (multi-trip)
    SOLVE_SECONDS     = 15


# ==============================================================================
# 2) DỮ LIỆU  (MẪU tổng hợp — thay bằng dữ liệu thật ở các hàm bên dưới)
# ==============================================================================
# Toạ độ gần Đà Nẵng: KS = bãi Khánh Sơn; depot theo quận; điểm tập kết ngẫu nhiên.
KS_LATLNG   = (16.036, 108.166)                    # Khánh Sơn (xấp xỉ)
DEPOTS      = {                                     # D = {depot theo quận}
    "Depot_HaiChau":  (16.052, 108.220),
    "Depot_LienChieu":(16.075, 108.150),
}

def sinh_diem_tap_ket(n=22, seed=7):
    """Sinh N điểm tập kết mẫu + lượng rác q_i (kg). THAY bằng dữ liệu thật/ML."""
    random.seed(seed)
    pts = []
    for i in range(n):
        lat = 16.036 + random.uniform(-0.05, 0.06)
        lng = 108.19 + random.uniform(-0.06, 0.05)
        qi  = random.choice([300, 500, 700, 900, 1100, 1400, 1800])  ### <<< q_i từ ML
        si  = random.choice([3, 4, 5, 6])                            # s_i (phút gom)
        pts.append({"name": f"n{i+1}", "latlng": (lat, lng), "q": qi, "s": si})
    return pts

def haversine_m(a, b):
    R = 6371000.0
    la1, lo1, la2, lo2 = map(math.radians, [a[0], a[1], b[0], b[1]])
    dla, dlo = la2 - la1, lo2 - lo1
    h = math.sin(dla/2)**2 + math.cos(la1)*math.cos(la2)*math.sin(dlo/2)**2
    return 2 * R * math.asin(math.sqrt(h))

# ==============================================================================
# 3) TIỀN XỬ LÝ (Hybrid): tách điểm ảo (dummy) nếu q_i > Q_max  (C3)
# ==============================================================================
def tien_xu_ly(diem, Qmax):
    """f(i) = ceil(q_i / Q_max). Nếu q_i > Q_max -> tách thành nhiều điểm ảo."""
    out = []
    for d in diem:
        f = max(1, math.ceil(d["q"] / Qmax))
        if f == 1:
            out.append(dict(d, q=d["q"], parent=d["name"]))
        else:
            chia = d["q"] / f
            for k in range(f):
                out.append(dict(d, name=f'{d["name"]}#{k+1}', q=chia, parent=d["name"]))
    return out


# ==============================================================================
# 4) DỰNG MÔ HÌNH OR-TOOLS
# ==============================================================================
def build_and_solve():
    Qmax = int(P.Q * P.THETA)                     # C1: sức chứa hiệu dụng θ·Q
    diem = tien_xu_ly(sinh_diem_tap_ket(), Qmax)  # danh sách điểm (đã tách dummy)
    C    = len(diem)

    # --- C10: số xe khả dụng = floor(K_max·(1 - spare)); phần còn lại là standby ---
    K_avail   = int(math.floor(P.K_MAX * (1 - P.SPARE_RATIO)))
    K_standby_reserved = P.K_MAX - K_avail
    V = K_avail

    # --- Gán mỗi xe cho 1 depot (xoay vòng giữa các depot theo quận) ---
    depot_names = list(DEPOTS.keys())
    veh_depot   = [depot_names[k % len(depot_names)] for k in range(V)]

    # -------- Danh sách NODE --------
    # [0..C-1]            : điểm tập kết (collection)
    # [C..C+R-1]          : bản sao KS trung gian (reload / multi-trip)  -> tuỳ chọn
    # [.. per-vehicle ..] : START copy (tại depot xe) + END copy (tại KS) cho từng xe
    R = P.RELOAD_COPIES
    node_loc, node_type, node_service, node_demand, node_label = [], [], [], [], []

    for d in diem:                                 # collection
        node_loc.append(d["latlng"]); node_type.append("N")
        node_service.append(d["s"]);  node_demand.append(int(round(d["q"])))
        node_label.append(d["name"])
    KS_reload_ids = []
    for r in range(R):                             # reload copies (KS trung gian)
        KS_reload_ids.append(len(node_loc))
        node_loc.append(KS_LATLNG); node_type.append("KS_RELOAD")
        node_service.append(P.S_KS + P.DT_KS); node_demand.append(-Qmax)
        node_label.append(f"KS_reload_{r+1}")
    start_ids, end_ids = [], []
    for k in range(V):                             # start (depot) + end (KS) mỗi xe
        start_ids.append(len(node_loc))
        node_loc.append(DEPOTS[veh_depot[k]]); node_type.append("START")
        node_service.append(0); node_demand.append(0)
        node_label.append(f"START[{veh_depot[k]}]")
    for k in range(V):
        end_ids.append(len(node_loc))
        node_loc.append(KS_LATLNG); node_type.append("END_KS")
        node_service.append(P.S_KS + P.DT_KS); node_demand.append(-Qmax)
        node_label.append("END_KS")
    Nn = len(node_loc)

    # -------- Ma trận khoảng cách d_ij (mét) & thời gian t_ij (phút) --------
    ### <<< CẮM DỮ LIỆU THẬT: thay 2 ma trận dưới bằng Google Distance Matrix / GPS
    dist = [[0]*Nn for _ in range(Nn)]
    tmin = [[0]*Nn for _ in range(Nn)]
    for i in range(Nn):
        for j in range(Nn):
            if i == j: continue
            dm = haversine_m(node_loc[i], node_loc[j]) * 1.3   # 1.3: hệ số đường bộ
            dist[i][j] = int(round(dm))
            tmin[i][j] = int(round(dm/1000.0 / P.SPEED_KMH * 60))

    # -------- OR-Tools manager & model --------
    mgr = pywrapcp.RoutingIndexManager(Nn, V, start_ids, end_ids)
    routing = pywrapcp.RoutingModel(mgr)

    # Arc cost = α·d_ij + β·(t_ij + s_i)
    #   -> khớp hàm mục tiêu: β·ΣT gồm cả Σt_ij, Σs_i và n_k^KS·(s_KS+Δt_KS).
    #   -> mỗi lần ghé KS tốn β·(s_KS+Δt_KS) nên KHÔNG còn ghé KS "rỗng" miễn phí.
    def arc_cost(fi, ti):
        i, j = mgr.IndexToNode(fi), mgr.IndexToNode(ti)
        return P.ALPHA * dist[i][j] + P.BETA * (tmin[i][j] + node_service[i])
    cost_cb = routing.RegisterTransitCallback(arc_cost)
    routing.SetArcCostEvaluatorOfAllVehicles(cost_cb)

    # γ·K_active : chi phí cố định mỗi xe được dùng
    routing.SetFixedCostOfAllVehicles(P.GAMMA)

    # ----- Dimension TẢI TRỌNG (C1, C5 mandatory dump reset) -----
    def demand_cb(fi):
        return node_demand[mgr.IndexToNode(fi)]
    dcb = routing.RegisterUnaryTransitCallback(demand_cb)
    routing.AddDimension(dcb, Qmax, Qmax, True, "Cap")   # slack=Qmax cho phép reset ở KS
    # (tải reset khi ghé KS nhờ demand âm -Qmax + cumul bị kẹp [0, Qmax])

    # ----- Dimension THỜI GIAN (C2 soft TW, C6 service time) -----
    def time_cb(fi, ti):
        i, j = mgr.IndexToNode(fi), mgr.IndexToNode(ti)
        return tmin[i][j] + node_service[i]              # + s_i khi rời i (C6)
    tcb = routing.RegisterTransitCallback(time_cb)
    routing.AddDimension(tcb, 120, P.HORIZON, True, "Time")   # fix_start=True -> đi lúc 20:00
    time_dim = routing.GetDimensionOrDie("Time")

    # C2: phạt mềm P·e_k nếu xe kết thúc (tại END_KS) sau 24:00 (=240')
    for k in range(V):
        end_idx = routing.End(k)
        time_dim.SetCumulVarSoftUpperBound(end_idx, P.SHIFT_END_MIN, P.PENALTY)

    # C3: reload copies là tuỳ chọn (được phép bỏ nếu không cần) — penalty 0
    for rid in KS_reload_ids:
        routing.AddDisjunction([mgr.NodeToIndex(rid)], 0)

    # ----- Cấu hình tìm kiếm -----
    prm = pywrapcp.DefaultRoutingSearchParameters()
    prm.first_solution_strategy = routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC
    prm.local_search_metaheuristic = routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH
    prm.time_limit.FromSeconds(P.SOLVE_SECONDS)

    sol = routing.SolveWithParameters(prm)
    ctx = dict(mgr=mgr, routing=routing, sol=sol, time_dim=time_dim,
               node_label=node_label, node_type=node_type, node_loc=node_loc,
               node_demand=node_demand, dist=dist, tmin=tmin, V=V, Qmax=Qmax,
               veh_depot=veh_depot, K_avail=K_avail,
               K_standby_reserved=K_standby_reserved)
    return ctx


# ==============================================================================
# 5) TRÍCH KẾT QUẢ + MODULE GATE BATCHING HẬU KỲ (C8, C9)
# ==============================================================================
def hhmm(t):
    """phút từ 20:00 -> 'HH:MM'."""
    total = 20*60 + int(round(t))
    return f"{(total//60)%24:02d}:{total%60:02d}"

def trich_ket_qua(ctx):
    mgr, routing, sol, time_dim = ctx["mgr"], ctx["routing"], ctx["sol"], ctx["time_dim"]
    cap_dim = routing.GetDimensionOrDie("Cap")
    routes, ks_events = [], []      # ks_events: (arrival_min, vehicle, node_label)
    for k in range(ctx["V"]):
        idx = routing.Start(k); seq = []
        used = False
        while not routing.IsEnd(idx):
            nd = mgr.IndexToNode(idx)
            arr = sol.Value(time_dim.CumulVar(idx))
            load = sol.Value(cap_dim.CumulVar(idx))
            seq.append((ctx["node_label"][nd], ctx["node_type"][nd], arr, load, nd))
            if ctx["node_type"][nd] == "N": used = True
            if ctx["node_type"][nd] == "KS_RELOAD":
                ks_events.append([arr, k, ctx["node_label"][nd]])
            idx = sol.Value(routing.NextVar(idx))
        # điểm END (KS chốt ca)
        nd = mgr.IndexToNode(idx); arr = sol.Value(time_dim.CumulVar(idx))
        seq.append((ctx["node_label"][nd], "END_KS", arr, 0, nd))
        if used: ks_events.append([arr, k, "END_KS"])
        routes.append({"veh": k, "depot": ctx["veh_depot"][k], "used": used, "seq": seq,
                       "end_min": arr})
    return routes, ks_events

def gate_batching(ks_events):
    """
    C8: xếp các lượt xe tới cổng KS thành NHÓM <= B_KS, mỗi nhóm mất B*t_gate,
        giữa 2 nhóm cách >= Δt_gap. Xe tới sớm chờ ngoài cổng -> w_k^KS.
    C9: nếu w_k^KS + (s_KS + Δt_KS) > T_timeout -> đánh dấu "SỰ CỐ".
    (Heuristic hậu kỳ: định tuyến đã xong, đây là lớp lập lịch cổng.)
    """
    ev = sorted(ks_events, key=lambda x: x[0])       # theo giờ tới cổng
    results, incidents = [], []
    i, bi, gate_free = 0, 0, None                    # gate 1 cửa, phục vụ theo đợt
    while i < len(ev):
        bi += 1
        first_arr = ev[i][0]
        # đợt bắt đầu khi cổng rảnh (>= đợt trước + Δt_gap) VÀ có xe đang chờ
        start = first_arr if gate_free is None else max(first_arr, gate_free + P.DT_GAP)
        grp = []                                     # gom tối đa B_KS xe đã tới trước 'start'
        j = i
        while j < len(ev) and len(grp) < P.B_KS and ev[j][0] <= start + 1e-9:
            grp.append(ev[j]); j += 1
        proc = len(grp) * P.T_GATE                   # thời gian qua cổng cả đợt
        gate_free = start + proc
        for (arr, k, lbl) in grp:
            w = start - arr                          # w_k^KS: chờ ngoài cổng
            tong = w + (P.S_KS + P.DT_KS)            # tổng chờ + lưu bãi
            row = dict(veh=k, node=lbl, arrive=arr, batch=bi, wait=w,
                       gate_start=start, leave=gate_free + (P.S_KS + P.DT_KS),
                       tong_cho_bai=tong, su_co=tong > P.T_TIMEOUT)
            results.append(row)
            if row["su_co"]: incidents.append(row)
        i = j
    return results, incidents


# ==============================================================================
# 6) IN KẾT QUẢ + VẼ BẢN ĐỒ
# ==============================================================================
def bao_cao(ctx, routes, ks_events):
    mgr, routing, sol = ctx["mgr"], ctx["routing"], ctx["sol"]
    print("="*78)
    print(" KẾT QUẢ TỐI ƯU — TUYẾN XE THU GOM RÁC (VYLT 2026)")
    print("="*78)
    if sol is None:
        print(" !! Không tìm được lời giải khả thi. Nới ΔT_max / tăng K_max / giảm q_i."); return

    used = [r for r in routes if r["used"]]
    K_active = len(used)
    tot_dist = tot_tmove = 0
    for k in range(ctx["V"]):
        idx = routing.Start(k)
        while not routing.IsEnd(idx):
            nx = sol.Value(routing.NextVar(idx))
            i, j = mgr.IndexToNode(idx), mgr.IndexToNode(nx)
            tot_dist += ctx["dist"][i][j]; tot_tmove += ctx["tmin"][i][j]
            idx = nx

    print(f" Số xe điều động  K_active = {K_active} / khả dụng {ctx['K_avail']} "
          f"(K_max={P.K_MAX}, standby dự phòng={ctx['K_standby_reserved']})")
    print(f" Tổng quãng đường Σd_ij   = {tot_dist/1000:.2f} km")
    print(f" Tổng t.g di chuyển Σt_ij = {tot_tmove} phút")
    print("-"*78)
    for r in used:
        print(f"\n  ▶ XE {r['veh']+1}  (xuất phát {r['depot']})")
        trip = 1; line = "     20:00 " + r["depot"]
        for (lbl, typ, arr, load, nd) in r["seq"]:
            if typ == "START":  continue
            if typ in ("KS_RELOAD", "END_KS"):
                line += f" → [{hhmm(arr)} KS ⭳đổ]"
                print(line)
                if typ == "KS_RELOAD":
                    trip += 1; line = f"     (chuyến {trip})"
                else:
                    line = ""
            else:
                line += f" → {hhmm(arr)} {lbl}({load}kg)"
        print(f"     ⇒ Kết thúc ca tại KS lúc {hhmm(r['end_min'])}"
              + ("  ⚠ TRỄ >24:00" if r['end_min'] > P.SHIFT_END_MIN else ""))

    # -------- Gate batching (C8/C9) --------
    gb, inc = gate_batching(ks_events)
    print("\n" + "-"*78)
    print(" LỊCH QUA CỔNG KS (C8 Gate Batching, ≤{} xe/nhóm, giãn {}′):".format(P.B_KS, P.DT_GAP))
    for row in sorted(gb, key=lambda x: (x["batch"], x["arrive"])):
        flag = "  ⚠ SỰ CỐ (>T_timeout)" if row["su_co"] else ""
        print(f"   Nhóm {row['batch']} | Xe {row['veh']+1:>2} | tới {hhmm(row['arrive'])}"
              f" | chờ w={row['wait']:.1f}′ | vào cổng {hhmm(row['gate_start'])}{flag}")
    print(f"\n C9: số lượt bị đánh dấu SỰ CỐ (chờ+bãi > {P.T_TIMEOUT}′) = {len(inc)}"
          + (" → cần xe standby định tuyến lại." if inc else " → không có."))
    print("="*78)
    return K_active, tot_dist, tot_tmove

def ve_ban_do(ctx, routes, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 8))
    cols = plt.cm.tab10.colors
    # KS + depot
    ax.scatter(*KS_LATLNG[::-1], c="red", marker="*", s=420, zorder=5, label="KS (Khánh Sơn)")
    for name,(la,ln) in DEPOTS.items():
        ax.scatter(ln, la, c="black", marker="s", s=130, zorder=5)
        ax.annotate(name, (ln, la), fontsize=8, xytext=(4,4), textcoords="offset points")
    for r in routes:
        if not r["used"]: continue
        c = cols[r["veh"] % 10]
        xs = [ctx["node_loc"][nd][1] for (_,_,_,_,nd) in r["seq"]]
        ys = [ctx["node_loc"][nd][0] for (_,_,_,_,nd) in r["seq"]]
        ax.plot(xs, ys, "-", color=c, lw=1.3, alpha=0.8, zorder=2, label=f"Xe {r['veh']+1}")
        for (lbl, typ, arr, load, nd) in r["seq"]:
            if typ == "N":
                ax.scatter(ctx["node_loc"][nd][1], ctx["node_loc"][nd][0],
                           color=c, s=26, zorder=3)
    ax.set_title("Tuyến thu gom rác tối ưu (MD-OVRPTW + VRPMD + Multi-trip)")
    ax.set_xlabel("Kinh độ"); ax.set_ylabel("Vĩ độ")
    ax.legend(fontsize=7, loc="best", ncol=2); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(path, dpi=130)
    print(f" [Đã lưu bản đồ tuyến: {path}]")


if __name__ == "__main__":
    ctx = build_and_solve()
    routes, ks_events = trich_ket_qua(ctx) if ctx["sol"] else ([], [])
    bao_cao(ctx, routes, ks_events)
    if ctx["sol"]:
        try:
            ve_ban_do(ctx, routes, "tuyen_thu_gom_rac.png")
        except Exception as e:
            print(" (Bỏ qua vẽ bản đồ:", e, ")")
