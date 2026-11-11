#!/usr/bin/env python3
"""P2-30 stage 4a (deliverable 10): one column list per log file. The header is built from the list and every row writer
counts its own fields against it. Each row keeps its own << chain and formats: `ofs <<` -> `row <<` on the writer's row,
assembled in a stream carrying the writer's format (copyfmt), counted, then written. Every replacement must match exactly
once or nothing is written."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
D = os.path.join(ROOT, "src", "Integrator", "MMWSpalling")

SHARED_INSERT_BEFORE = """    // =====================================================================
    // P2-18 A3: THE STEP PROLOGUE both paths share (formerly the first lines
"""
SHARED_NEW = """    // =====================================================================
    // P2-30: ONE COLUMN LIST PER LOG FILE. The ordered column NAMES of each CSV this class writes live here; the header is
    // built from the list, and every row writer counts its own fields against it before writing (it assembles the row in
    // a local stream carrying the writer's own format, counts, aborts on a mismatch, then writes). A row writer keeps its
    // own << chain and formats (a column can print as `0` in one row type and `0.0` in another). This catches a missing
    // or an extra field, NOT two fields swapped.
    // =====================================================================
    static std::string CsvHeader(const std::vector<std::string>& cols)
    {
        std::string s;
        for (std::size_t i = 0; i < cols.size(); ++i)
        {
            if (i) s += ",";
            s += cols[i];
        }
        return s + "\\n";
    }

    // `row` is one finished row (its trailing newline is not a field); every field is a number, so a comma is a separator
    static void CsvCheckRow(const std::string& row, std::size_t n_cols, const char* file)
    {
        std::size_t n = 1;
        for (char ch : row) if (ch == ',') ++n;
        if (n != n_cols)
            Util::Abort(INFO, file, ": a row has ", n, " fields, its header has ", n_cols, ": ", row);
    }

    // <plot_file>_removal_events.csv: every row writer under the face list (firing, vapour, feet) shares this one list
    static std::vector<std::string> FrEventsColumns(bool bk_on)
    {
        std::vector<std::string> c = {"time", "col_i", "col_j", "regime", "k_top", "T_top", "a_f", "Sp_top",
                                      "h_scan", "h_applied", "n_voided",
                                      "T_surf", "ell_eff", "h_f", "T_bar_flake",
                                      "rc_flag", "d_ell_rel", "ell_raw", "delta", "q_ratio",
                                      // P2-9 deliverable 6, appended: the line
                                      "dir", "mu", "k", "run", "h_line", "h_n", "f_after", "score", "ts_fb", "sheet",
                                      "through", "q_net", "k_c"};
        if (bk_on) c.push_back("bk");                      // P2-15: under spall.buckle_span > 0 only, the line's ratio
        // P2-19, appended: the line's profile depths and the energy thin limit [m] (diagnostics; nothing reads them)
        for (const char* s : {"ell_p", "z_peak", "ell_pk", "h_E"}) c.push_back(s);
        return c;
    }

    // <plot_file>_buckle_log.csv
    static std::vector<std::string> FrBuckleColumns()
    {
        return {"time", "col_i", "col_j", "class", "fire", "through", "fb", "nodepth", "score", "h_col",
                "T_surf", "ell_p", "z_peak", "ell_pk", "ell_eff", "x", "bk", "bk_lb", "blocked"};
    }

    // <plot_file>_jet_profile.csv (the optional columns follow the march's own switches)
    static std::vector<std::string> JetProfileColumns(bool expo, bool fsurf, bool count)
    {
        std::vector<std::string> c = {"time", "bin", "r_lo", "n_cols", "T_gas", "P_bin", "P_share"};
        if (expo) c.push_back("T_eff");
        if (fsurf) { c.push_back("n_fs"); c.push_back("m_ratio"); }
        if (count) c.push_back("P_neg");
        return c;
    }

    // <plot_file>_side_power.csv
    static std::vector<std::string> SidePowerColumns()
    {
        return {"time", "col_i", "col_j", "r", "P_side"};
    }

"""

EDITS = {
    "Shared.H": [(SHARED_INSERT_BEFORE, SHARED_NEW + SHARED_INSERT_BEFORE)],
    "FaceRemoval.H": [
        # the event log's one header: built from the list
        ('''    static std::string FrEventsHeader(bool bk_on)
    {
        std::ostringstream ofs;
        ofs << "time,col_i,col_j,regime,k_top,T_top,a_f,Sp_top,"
               "h_scan,h_applied,n_voided"
            << ",T_surf,ell_eff,h_f,T_bar_flake,"
               "rc_flag,d_ell_rel,ell_raw,delta,q_ratio"
            // P2-9 deliverable 6, appended: the line
            << ",dir,mu,k,run,h_line,h_n,f_after,score,ts_fb,sheet,"
               "through,q_net,k_c"
            // P2-15: under spall.buckle_span > 0 only, the line's ratio
            << (bk_on ? ",bk" : "")
            // P2-19, appended: the line's profile depths and the energy
            // thin limit [m] (diagnostics; nothing reads them)
            << ",ell_p,z_peak,ell_pk,h_E"
            << "\\n";
        return ofs.str();
    }
''', '''    static std::string FrEventsHeader(bool bk_on)
    {
        return CsvHeader(FrEventsColumns(bk_on));
    }
'''),
        # the feet row
        ('''                const Set::Scalar nan = std::numeric_limits<Set::Scalar>::quiet_NaN();
                for (int idx = 0; idx < ncols; ++idx)
                {
                    if (n_cut[idx] == 0) continue;
                    const Set::Scalar h = Set::Scalar(n_cut[idx]) * dz;
                    const int kt = k_top_global[static_cast<std::size_t>(idx)];
                    ofs << time << "," << (idx % nx_dom + dom_ilo) << "," << (idx / nx_dom + dom_jlo) << ",4,"
                        << kt << ",0,0,0,0," << h << "," << flip[idx]
                        << ",0,0,0,0,0,0,0,0,0";
                    ofs << ",5,1," << kt << ",0," << h << "," << h << ",0,0,0,0,0,0,0";
                    if (bk_on) ofs << ",0";
                    ofs << "," << nan << "," << nan << "," << nan << "," << nan << "\\n";
                }
''', '''                const Set::Scalar nan = std::numeric_limits<Set::Scalar>::quiet_NaN();
                const std::size_t n_ev = FrEventsColumns(bk_on).size();
                for (int idx = 0; idx < ncols; ++idx)
                {
                    if (n_cut[idx] == 0) continue;
                    const Set::Scalar h = Set::Scalar(n_cut[idx]) * dz;
                    const int kt = k_top_global[static_cast<std::size_t>(idx)];
                    std::ostringstream row;
                    row.copyfmt(ofs);
                    row << time << "," << (idx % nx_dom + dom_ilo) << "," << (idx / nx_dom + dom_jlo) << ",4,"
                        << kt << ",0,0,0,0," << h << "," << flip[idx]
                        << ",0,0,0,0,0,0,0,0,0";
                    row << ",5,1," << kt << ",0," << h << "," << h << ",0,0,0,0,0,0,0";
                    if (bk_on) row << ",0";
                    row << "," << nan << "," << nan << "," << nan << "," << nan << "\\n";
                    CsvCheckRow(row.str(), n_ev, "_removal_events.csv");
                    ofs << row.str();
                }
'''),
        # the buckle log's header and row
        ('''                    if (need_header)
                        ofs << "time,col_i,col_j,class,fire,through,fb,nodepth,score,h_col,"
                               "T_surf,ell_p,z_peak,ell_pk,ell_eff,x,bk,bk_lb,blocked\\n";
                    ofs.setf(std::ios::scientific);
                    ofs.precision(16);
                    for (const BkRow& br : bk_rows)
                    {
                        const FrLine& L = lines[br.l];
                        const LineOut& o = out[br.l];
                        const FrCell& c = cells[static_cast<std::size_t>(L.cell)];
                        ofs << time << "," << c.i << "," << c.j << "," << FrClass(L.dir)
                            << ",1," << br.through << "," << o.fb << "," << o.nodepth
                            << "," << o.score << "," << o.h_col << "," << o.rec.T_surf
                            << "," << o.ell_p << "," << o.z_peak << "," << o.ell_pk
                            << "," << o.rec.ell_eff << "," << br.x << "," << br.bk
                            << "," << br.bk_lb << "," << br.blocked << "\\n";
                    }
''', '''                    if (need_header) ofs << CsvHeader(FrBuckleColumns());
                    ofs.setf(std::ios::scientific);
                    ofs.precision(16);
                    const std::size_t n_bk = FrBuckleColumns().size();
                    for (const BkRow& br : bk_rows)
                    {
                        const FrLine& L = lines[br.l];
                        const LineOut& o = out[br.l];
                        const FrCell& c = cells[static_cast<std::size_t>(L.cell)];
                        std::ostringstream row;
                        row.copyfmt(ofs);
                        row << time << "," << c.i << "," << c.j << "," << FrClass(L.dir)
                            << ",1," << br.through << "," << o.fb << "," << o.nodepth
                            << "," << o.score << "," << o.h_col << "," << o.rec.T_surf
                            << "," << o.ell_p << "," << o.z_peak << "," << o.ell_pk
                            << "," << o.rec.ell_eff << "," << br.x << "," << br.bk
                            << "," << br.bk_lb << "," << br.blocked << "\\n";
                        CsvCheckRow(row.str(), n_bk, "_buckle_log.csv");
                        ofs << row.str();
                    }
'''),
        # the firing row
        ('''                    if (need_header) ofs << FrEventsHeader(bk_on);
                    ofs.setf(std::ios::scientific);
                    ofs.precision(16);
                    for (const Ev& ev : evs)
                    {
''', '''                    if (need_header) ofs << FrEventsHeader(bk_on);
                    ofs.setf(std::ios::scientific);
                    ofs.precision(16);
                    const std::size_t n_ev = FrEventsColumns(bk_on).size();
                    for (const Ev& ev : evs)
                    {
'''),
        ('''                        ofs << time << "," << c.i << "," << c.j << ",1,"
                            << k_top_global[static_cast<std::size_t>(c.col)] << ","
                            << Gl[0] << "," << Gc[2 * ev.l] << "," << 0.0 << ","
                            << o.h_col << "," << (eval_only ? 0.0 : ev.h_line) << "," << ev.n_v
                            << "," << o.rec.T_surf << ","
                            << (std::isfinite(o.rec.ell_eff) ? o.rec.ell_eff : 0.0)
                            << "," << o.rec.h_f << "," << o.rec.T_bar << ","
                            << Set::Scalar(flag) << ",0,"
                            << (std::isfinite(o.rec.ell_eff) ? o.rec.ell_eff : 0.0)
                            << ",0,0";
                        ofs << "," << L.dir << "," << L.mu << "," << c.k << "," << o.run
                            << "," << ev.h_line << "," << L.mu * ev.h_line
                            << "," << ev.f_after << "," << o.score << "," << o.fb
                            << "," << c.sheet << "," << o.through << "," << o.q_net
                            << "," << Gc[2 * ev.l + 1];
                        if (bk_on) ofs << "," << o.bk;
                        ofs << "," << o.ell_p << "," << o.z_peak << "," << o.ell_pk
                            << "," << o.h_E;
                        ofs << "\\n";
                    }
''', '''                        std::ostringstream row;
                        row.copyfmt(ofs);
                        row << time << "," << c.i << "," << c.j << ",1,"
                            << k_top_global[static_cast<std::size_t>(c.col)] << ","
                            << Gl[0] << "," << Gc[2 * ev.l] << "," << 0.0 << ","
                            << o.h_col << "," << (eval_only ? 0.0 : ev.h_line) << "," << ev.n_v
                            << "," << o.rec.T_surf << ","
                            << (std::isfinite(o.rec.ell_eff) ? o.rec.ell_eff : 0.0)
                            << "," << o.rec.h_f << "," << o.rec.T_bar << ","
                            << Set::Scalar(flag) << ",0,"
                            << (std::isfinite(o.rec.ell_eff) ? o.rec.ell_eff : 0.0)
                            << ",0,0";
                        row << "," << L.dir << "," << L.mu << "," << c.k << "," << o.run
                            << "," << ev.h_line << "," << L.mu * ev.h_line
                            << "," << ev.f_after << "," << o.score << "," << o.fb
                            << "," << c.sheet << "," << o.through << "," << o.q_net
                            << "," << Gc[2 * ev.l + 1];
                        if (bk_on) row << "," << o.bk;
                        row << "," << o.ell_p << "," << o.z_peak << "," << o.ell_pk
                            << "," << o.h_E;
                        row << "\\n";
                        CsvCheckRow(row.str(), n_ev, "_removal_events.csv");
                        ofs << row.str();
                    }
'''),
        # the vapour row
        ('''                        const Set::Scalar h_v = vap_f[ci] * dx;
                        ofs << time << "," << c.i << "," << c.j << ",2,"
                            << k_top_global[static_cast<std::size_t>(c.col)] << ","
                            << vap_T[ci] << ",0,0,0," << h_v << ",1"
                            << ",0,0,0,0,0,0,0,0,0";
                        ofs << "," << fl_faces[static_cast<std::size_t>(c.qline)].dir << ",1," << c.k << ",0,"
                            << h_v << "," << h_v << ",0,0,0," << c.sheet << ",0,0,0";
                        if (bk_on) ofs << ",0";
                        ofs << "," << nan << "," << nan << "," << nan << "," << nan << "\\n";
''', '''                        const Set::Scalar h_v = vap_f[ci] * dx;
                        std::ostringstream row;
                        row.copyfmt(ofs);
                        row << time << "," << c.i << "," << c.j << ",2,"
                            << k_top_global[static_cast<std::size_t>(c.col)] << ","
                            << vap_T[ci] << ",0,0,0," << h_v << ",1"
                            << ",0,0,0,0,0,0,0,0,0";
                        row << "," << fl_faces[static_cast<std::size_t>(c.qline)].dir << ",1," << c.k << ",0,"
                            << h_v << "," << h_v << ",0,0,0," << c.sheet << ",0,0,0";
                        if (bk_on) row << ",0";
                        row << "," << nan << "," << nan << "," << nan << "," << nan << "\\n";
                        CsvCheckRow(row.str(), n_ev, "_removal_events.csv");
                        ofs << row.str();
'''),
    ],
    "Flame.H": [
        ('''            if (prof)
            {
                prow << time << "," << ib << "," << Set::Scalar(ib) * dr << "," << bin.size()
                     << "," << T_in << "," << P_bin << "," << share;
                if (expo)  prow << "," << T_col;
                if (fsurf) prow << "," << n_fs << "," << m_cur / m_jet;
                if (count) prow << "," << P_nb * dA;
                prow << "\\n";
            }
''', '''            if (prof)
            {
                std::ostringstream row;
                row.copyfmt(prow);
                row << time << "," << ib << "," << Set::Scalar(ib) * dr << "," << bin.size()
                     << "," << T_in << "," << P_bin << "," << share;
                if (expo)  row << "," << T_col;
                if (fsurf) row << "," << n_fs << "," << m_cur / m_jet;
                if (count) row << "," << P_nb * dA;
                row << "\\n";
                CsvCheckRow(row.str(), JetProfileColumns(expo, fsurf, count).size(), "_jet_profile.csv");
                prow << row.str();
            }
'''),
        ('''                if (first)
                {
                    ofs << "time,bin,r_lo,n_cols,T_gas,P_bin,P_share";
                    if (expo)  ofs << ",T_eff";
                    if (fsurf) ofs << ",n_fs,m_ratio";
                    if (count) ofs << ",P_neg";
                    ofs << "\\n";
                }
''', '''                if (first) ofs << CsvHeader(JetProfileColumns(expo, fsurf, count));
'''),
    ],
    "Ledger.H": [
        ('''                            if (first) ofs << "time,col_i,col_j,r,P_side\\n";
                            ofs.setf(std::ios::scientific);
                            ofs.precision(10);
''', '''                            if (first) ofs << CsvHeader(SidePowerColumns());
                            ofs.setf(std::ios::scientific);
                            ofs.precision(10);
                            const std::size_t n_side = SidePowerColumns().size();
'''),
        ('''                                ofs << time << "," << (ci + ilo_dom) << ","
                                    << (cj + jlo_dom) << ","
                                    << std::sqrt(cx * cx + cy * cy) << ","
                                    << sidecol[c] << "\\n";
''', '''                                std::ostringstream row;
                                row.copyfmt(ofs);
                                row << time << "," << (ci + ilo_dom) << ","
                                    << (cj + jlo_dom) << ","
                                    << std::sqrt(cx * cx + cy * cy) << ","
                                    << sidecol[c] << "\\n";
                                CsvCheckRow(row.str(), n_side, "_side_power.csv");
                                ofs << row.str();
'''),
    ],
}


def main():
    new = {}
    for fn, reps in EDITS.items():
        p = os.path.join(D, fn)
        t = open(p).read()
        for old, rep in reps:
            n = t.count(old)
            if n != 1:
                sys.exit("%s: expected 1 match, found %d of:\n%s" % (fn, n, old[:220]))
            t = t.replace(old, rep, 1)
        new[p] = t
    for p, t in new.items():
        open(p, "w").write(t)
    print("4a applied to", ", ".join(sorted(EDITS)))


if __name__ == "__main__":
    main()
