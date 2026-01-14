import gurobipy as gp
from gurobipy import GRB
import networkx as nx
import heapq
from dataclasses import dataclass
from typing import Dict, Tuple, List, Optional

@dataclass(frozen=True)
class Demand:
    o: int
    d: int
    q: float
    Tmax: float
    taumax: int

@dataclass(frozen=True)
class Route:
    nodes: Tuple[int, ...]
    edges: Tuple[Tuple[int, int], ...]
    time: float
    transfers: int
    is_direct: bool

class CG_LPP_SL:
    def __init__(
        self,
        G: nx.DiGraph,
        demands: Dict[str, Demand],
        alpha_direct: float,
        beta_direct_share: float,
        cap_per_freq: float,
        cost_per_freq_edge: float,
        fleet_limit: float,
        w_time: float,
        w_transfer: float,
        direct_max_transfers: int = 1,
        outputflag: int = 0
    ):
        self.G = G
        self.demands = demands
        self.alpha_direct = alpha_direct
        self.beta_direct_share = beta_direct_share
        self.cap_per_freq = cap_per_freq
        self.cost_per_freq_edge = cost_per_freq_edge
        self.fleet_limit = fleet_limit
        self.w_time = w_time
        self.w_transfer = w_transfer
        self.direct_max_transfers = direct_max_transfers
        self.outputflag = outputflag

        for u, v in self.G.edges():
            if "time" not in self.G[u][v]:
                raise ValueError("Fiecare arc trebuie să aibă atributul 'time' (numeric).")

        self.L = list(self.G.edges())
        self.edge_time = {(u, v): float(self.G[u][v]["time"]) for (u, v) in self.L}

        self.shortest_time = {}
        for k, dem in self.demands.items():
            self.shortest_time[k] = nx.shortest_path_length(self.G, dem.o, dem.d, weight="time")

        self.P: Dict[str, List[Route]] = {k: [] for k in self.demands}

    def _build_route(self, k: str, nodes: List[int]) -> Route:
        dem = self.demands[k]
        edges = [(nodes[i], nodes[i + 1]) for i in range(len(nodes) - 1)]
        t = sum(self.edge_time[e] for e in edges)
        transfers = max(0, len(edges) - 1)
        is_direct = (transfers <= self.direct_max_transfers) and (t <= self.alpha_direct * self.shortest_time[k])
        return Route(nodes=tuple(nodes), edges=tuple(edges), time=t, transfers=transfers, is_direct=is_direct)

    def init_pool_with_shortest_paths(self):
        for k, dem in self.demands.items():
            nodes = nx.shortest_path(self.G, dem.o, dem.d, weight="time")
            r = self._build_route(k, nodes)
            if r.time > dem.Tmax or r.transfers > dem.taumax:
                raise RuntimeError(
                    f"Ruta inițială shortest-path pentru {k} încalcă service-level. "
                    f"Crește Tmax/taumax sau modifică graful."
                )
            self.P[k].append(r)

    def solve_rmp_lp(self):
        m = gp.Model("RMP_LP")
        m.setParam("OutputFlag", self.outputflag)

        x = {}
        for (u, v) in self.L:
            x[(u, v)] = m.addVar(lb=0.0, vtype=GRB.CONTINUOUS, obj=self.cost_per_freq_edge, name=f"x_{u}_{v}")

        y = {}
        for k in self.P:
            dem = self.demands[k]
            for idx, r in enumerate(self.P[k]):
                cost_route_unit = self.w_time * r.time + self.w_transfer * r.transfers
                y[(k, idx)] = m.addVar(
                    lb=0.0, ub=1.0, vtype=GRB.CONTINUOUS,
                    obj=dem.q * cost_route_unit,
                    name=f"y_{k}_{idx}"
                )

        m.update()

        con_demand = {}
        for k in self.P:
            con_demand[k] = m.addConstr(
                gp.quicksum(y[(k, idx)] for idx in range(len(self.P[k]))) == 1.0,
                name=f"Demand_{k}"
            )

        con_cap = {}
        for e in self.L:
            expr = gp.LinExpr()
            for k in self.P:
                dem = self.demands[k]
                for idx, r in enumerate(self.P[k]):
                    if e in r.edges:
                        expr.addTerms(dem.q, y[(k, idx)])
            con_cap[e] = m.addConstr(expr <= self.cap_per_freq * x[e], name=f"Cap_{e[0]}_{e[1]}")

        con_fleet = m.addConstr(
            gp.quicksum(self.edge_time[e] * x[e] for e in self.L) <= self.fleet_limit,
            name="Fleet"
        )

        total_demand = sum(self.demands[k].q for k in self.demands)
        expr_direct = gp.LinExpr()
        for k in self.P:
            dem = self.demands[k]
            for idx, r in enumerate(self.P[k]):
                if r.is_direct:
                    expr_direct.addTerms(dem.q, y[(k, idx)])

        con_beta = m.addConstr(expr_direct >= self.beta_direct_share * total_demand, name="BetaDirect")

        m.optimize()
        if m.status != GRB.OPTIMAL:
            if m.status == GRB.INFEASIBLE:
                m.computeIIS()
                m.write("rmp_iis.ilp")
                print("RMP infeasible. IIS written to rmp_iis.ilp")
            raise RuntimeError(f"RMP LP not optimal. Status={m.status}")

        pi = {k: con_demand[k].Pi for k in con_demand}
        sigma = {e: con_cap[e].Pi for e in con_cap}
        top = sorted(sigma.items(), key=lambda kv: kv[1], reverse=True)[:5]
        print("Top sigma (capacity duals):", [(e, round(val, 4)) for e, val in top])

        gamma = con_fleet.Pi
        delta = con_beta.Pi

        return float(m.ObjVal), pi, sigma, gamma, delta

    def pricing_one_demand(
        self,
        k: str,
        pi_k: float,
        sigma: Dict[Tuple[int, int], float],
        delta: float,
        eps: float = 1e-7
    ) -> Optional[Route]:
        dem = self.demands[k]
        o, d = dem.o, dem.d

        Tmax = dem.Tmax
        taumax = dem.taumax
        Tshort = self.shortest_time[k]

        pq = []
        heapq.heappush(pq, (0.0, 0.0, 0, o, (o,)))
        best = {}

        best_route = None
        best_rc = 0.0

        while pq:
            rc_partial, time_used, edges_used, u, path_nodes = heapq.heappop(pq)
            transfers = max(0, edges_used - 1)
            if time_used > Tmax or transfers > taumax:
                continue

            key = (u, edges_used)
            if key in best and best[key] <= rc_partial:
                continue
            best[key] = rc_partial

            if u == d and edges_used >= 1:
                is_direct = (transfers <= self.direct_max_transfers) and (time_used <= self.alpha_direct * Tshort)

                c_pass = dem.q * (self.w_time * time_used + self.w_transfer * transfers)
                dual_cap = dem.q * sum(sigma[(path_nodes[i], path_nodes[i + 1])] for i in range(len(path_nodes) - 1))
                dual_beta = dem.q * (delta if is_direct else 0.0)

                rc = c_pass - pi_k + dual_cap - dual_beta

                if rc < best_rc - eps:
                    best_rc = rc
                    edges = tuple((path_nodes[i], path_nodes[i + 1]) for i in range(len(path_nodes) - 1))
                    best_route = Route(nodes=tuple(path_nodes), edges=edges, time=time_used, transfers=transfers, is_direct=is_direct)
                continue

            for v in self.G.successors(u):
                if v in path_nodes:
                    continue
                e = (u, v)

                new_time = time_used + self.edge_time[e]
                new_edges = edges_used + 1
                new_transfers = max(0, new_edges - 1)
                if new_time > Tmax or new_transfers > taumax:
                    continue

                arc_cost = dem.q * (self.w_time * self.edge_time[e]) + dem.q * sigma[e]
                new_rc_partial = rc_partial + arc_cost

                heapq.heappush(pq, (new_rc_partial, new_time, new_edges, v, path_nodes + (v,)))

        return best_route

    def run_column_generation(self, max_iter: int = 50, eps: float = 1e-7):
        self.init_pool_with_shortest_paths()

        total_demand = sum(self.demands[k].q for k in self.demands)
        direct_possible = 0.0
        for k, dem in self.demands.items():
            r0 = self.P[k][0]
            if r0.is_direct:
                direct_possible += dem.q
        print(f"Initial direct-share possible (from init pool): {direct_possible/total_demand:.3f} (target beta={self.beta_direct_share:.3f})")

        for it in range(1, max_iter + 1):
            obj, pi, sigma, gamma, delta = self.solve_rmp_lp()

            added = 0
            for k in self.demands:
                new_r = self.pricing_one_demand(k, pi[k], sigma, delta, eps=eps)
                if new_r is None:
                    continue
                if any(r.nodes == new_r.nodes for r in self.P[k]):
                    continue

                dem = self.demands[k]
                c_pass = dem.q * (self.w_time * new_r.time + self.w_transfer * new_r.transfers)
                dual_cap = dem.q * sum(sigma[e] for e in new_r.edges)
                dual_beta = dem.q * (delta if new_r.is_direct else 0.0)
                rc = c_pass - pi[k] + dual_cap - dual_beta

                if rc < -1e-6:
                    self.P[k].append(new_r)
                    added += 1

            print(f"Iter {it:02d} | RMP obj={obj:.4f} | added_routes={added} | delta(Beta)={delta:.6f} | gamma(Fleet)={gamma:.6f}")

            if added == 0:
                print("Convergență: nu există rute cu cost redus negativ.")
                break

def demo():
    G = nx.DiGraph()

    # Două coridoare 0->4:
    # Fast: 0-1-3-4 (timp mic)
    # Slow: 0-2-5-4 (timp mai mare)
    edges_data = [
        (0, 1, 2), (1, 3, 2), (3, 4, 2),      # FAST corridor
        (0, 2, 4), (2, 5, 4), (5, 4, 4),      # SLOW corridor

        # muchii de legătură (permit rute mixte / ocoliri)
        (1, 2, 1), (3, 5, 1),
        (2, 3, 3), (1, 5, 3),

        # cereri locale (ca să încarce arcele fast)
        (1, 4, 3), (0, 3, 3), (2, 4, 5),
    ]

    for u, v, t in edges_data:
        G.add_edge(u, v, time=t)
        G.add_edge(v, u, time=t)

    demands = {
        # împing puternic pe coridorul FAST (0-1-3-4)
        "k1": Demand(o=0, d=4, q=240, Tmax=20, taumax=4),
        "k2": Demand(o=1, d=4, q=180, Tmax=18, taumax=4),
        "k3": Demand(o=0, d=3, q=160, Tmax=15, taumax=4),

        # cereri care pot ocoli prin SLOW sau mixte
        "k4": Demand(o=2, d=4, q=120, Tmax=22, taumax=4),
        "k5": Demand(o=0, d=5, q=100, Tmax=22, taumax=4),
    }

    solver = CG_LPP_SL(
        G=G,
        demands=demands,
        alpha_direct=1.4,
        beta_direct_share=0.6,
        cap_per_freq=29.0,
        cost_per_freq_edge=12.0,
        fleet_limit=110.0,
        w_time=1.0,
        w_transfer=6.0,
        direct_max_transfers=2,
        outputflag=0
    )

    solver.run_column_generation(max_iter=30)


if __name__ == "__main__":
    demo()
