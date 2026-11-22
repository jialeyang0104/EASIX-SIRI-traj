"""The current Figures S6 and S7, built only from released aggregates."""
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import NullLocator, PercentFormatter
import numpy as np
import pandas as pd

COLORS = {1: "#4C78A8", 2: "#59A14F", 3: "#E15759"}
DB = {"mimiciv": "MIMIC-IV", "mimiciii": "MIMIC-III"}
ENDS = {"28-day": ("28-day mortality", "HR"), "ICU": ("ICU mortality", "sHR")}
LABELS = {
    "nee_any_24h": "Any NEE",
    "nee_dose": "NEE >0.30 vs 0<NEE≤0.30 µg/kg/min",
    "net_fluid": "Net fluid balance >35 vs ≤35 mL/kg",
    "documented_crrt_24h": "CRRT",
}
GRAY = "#53636E"


def ptext(p):
    if pd.isna(p):
        return "—"
    return "<0.001" if p < .001 else f"{p:.3f}"


def fixed_class_forest(aggregate_dir):
    forest_data = pd.read_csv(aggregate_dir / "figure_s6_fixed.csv")
    assert len(forest_data) == 48
    assert not forest_data.duplicated(["dataset", "outcome", "exposure", "trajectory_class"]).any()
    w,h=180,190
    fig=plt.figure(figsize=(w/25.4,h/25.4))
    def text(x,y,label,**kw):
        return fig.text(x/w,y/h,label,va="center",**kw)
    for row,(database,name) in enumerate(DB.items()):
        for col,(outcome,(endpoint,metric)) in enumerate(ENDS.items()):
            x,y0=3+90*col,182-92*row
            text(x,y0,chr(65+2*row+col),weight="bold",fontsize=10)
            text(x+6,y0,f"{name} · {endpoint}",weight="bold",fontsize=9.4)
            for i,(exposure,label) in enumerate(LABELS.items()):
                y=y0-9-18*i
                q=forest_data.query("dataset==@database and outcome==@outcome and exposure==@exposure").set_index("trajectory_class")
                assert len(q)==3 and q.interaction_p.nunique()==1
                text(x,y,label,weight="bold",fontsize=8.3)
                text(x+84,y,"P int "+ptext(q.interaction_p.iloc[0]),ha="right",color=GRAY,fontsize=8.3)
                ax=fig.add_axes([(x+12)/w,(y-11.7)/h,25/w,8.5/h])
                ax.set(xscale="log",xlim=(.2,6),ylim=(-.36,2.36),yticks=[])
                ax.axvline(1,color=GRAY,ls=":",lw=.8)
                ax.xaxis.set_minor_locator(NullLocator())
                ax.set_xticks([.2,1,6],labels=["0.2","1","6"])
                ax.tick_params(axis="x",length=2,pad=1,labelsize=8)
                ax.spines["left"].set_visible(False)
                if i<3:
                    ax.tick_params(axis="x",labelbottom=False,length=0)
                for k,marker in [(1,"o"),(2,"s"),(3,"^")]:
                    yy=3-k; r=q.loc[k]; c=COLORS[k]
                    ytext=y-11.7+(yy+.36)*8.5/2.72
                    text(x+10,ytext,f"C{k}",ha="right",color=c,fontsize=8.1)
                    if r.status!="estimated":
                        text(x+40,ytext,"Not estimable",color=GRAY,fontsize=8.1)
                        continue
                    ax.plot([max(.2,r.ci_low),min(6,r.ci_high)],[yy,yy],c=c,lw=1)
                    ax.plot(np.clip(r.ratio,.2,6),yy,marker=marker,c=c,ms=3.4,ls="none")
                    for flag,pos,m in [(r.ci_low<.2,.2,"<"),(r.ci_high>6,6,">")]:
                        if flag: ax.plot(pos,yy,c=c,marker=m,ms=3,clip_on=False)
                    text(x+40,ytext,f"{r.ratio:.2f} ({r.ci_low:.2f}–{r.ci_high:.2f})",fontsize=8.1)
                    text(x+84,ytext,ptext(r.p_value),ha="right",fontsize=8.1)
    text(w/2,5,"28-day mortality: HR · ICU mortality: sHR (95% CI); row P · P int: interaction P",
         ha="center",color=GRAY,fontsize=8.1)
    return fig


def treatment_probability(aggregate_dir):
    KEY = ["exposure", "dataset", "outcome"]
    curves = pd.read_csv(aggregate_dir / "figure_s7_four_treatment_curves.csv")
    tests = pd.read_csv(aggregate_dir / "figure_s7_four_treatment_tests.csv").set_index(KEY, verify_integrity=True)
    support = pd.read_csv(aggregate_dir / "figure_s7_four_treatment_histograms.csv")
    panels = [
        ("exposure_any_nee", "Any NEE"),
        ("exposure_nee_dose", "NEE dose"),
        ("exposure_net_fluid", "Net fluid balance"),
        ("exposure_crrt", "CRRT"),
    ]
    cohorts = [("mimiciv", "MIMIC-IV"), ("mimiciii", "MIMIC-III")]
    outcomes = [("28-day", "28-day mortality", "HR"), ("ICU", "ICU mortality", "sHR")]
    palette = {2: "#59A14F", 3: "#E15759"}
    line_styles = {2: "-", 3: (0, (3, 1.5))}
    gray = "#465864"
    light = "#DBE4E8"
    assert len(tests) == 16 and set(tests.status) <= {"estimated", "not_estimable"}
    assert set(tests.index) == {(t, d, e) for t, *_ in panels for d, _ in cohorts for e, *_ in outcomes}
    assert not curves.duplicated(KEY + ["predicted_class", "probability"]).any()
    assert curves.groupby(KEY + ["predicted_class"]).size().eq(101).all()
    assert (0 < curves.ci_low).all() and (curves.ci_low <= curves.ratio).all()
    assert (curves.ratio <= curves.ci_high).all() and np.isfinite(curves.ci_high).all()
    check = support.groupby(KEY + ["predicted_class", "exposure_arm"]).agg(
        n=("count", "sum"), arm_n=("arm_n", "first"), bins=("count", "size"))
    assert len(check) == 16 * 2 * 2 and check.bins.eq(20).all() and check.n.eq(check.arm_n).all()
    max_share = (support["count"] / support.arm_n.where(support.arm_n > 0)).max()

    width_mm, height_mm = 180, 203
    W, H = width_mm, height_mm
    fig = plt.figure(figsize=(width_mm / 25.4, height_mm / 25.4))


    def text(x, y, label, **kwargs):
        return fig.text(x / W, y / H, label, va="center", **kwargs)


    def axes_at(x, y, w, h):
        return fig.add_axes([x / W, y / H, w / W, h / H])


    def p_format(p):
        return "<0.001" if p < .001 else f"= {p:.3f}"


    def support_strip(x, y, w, rows, show_ticks):
        """For each class, show the lower/no arm above the higher/yes arm."""
        ax = axes_at(x, y, w, 4.3)
        patches, colors = [], []
        for r in rows.itertuples():
            lane = 3 - (2 * (r.predicted_class - 2) + r.exposure_arm)
            patches.append(Rectangle((r.bin_lower, lane + .10), r.bin_upper - r.bin_lower, .80))
            share = r.count / r.arm_n if r.arm_n else 0
            colors.append(to_rgba(palette[r.predicted_class], .05 + .9 * share / max_share))
        ax.add_collection(PatchCollection(patches, facecolors=colors, edgecolors="none"))
        ax.set(xlim=(0, 1), ylim=(0, 4), yticks=[], xticks=[0,.5,1] if show_ticks else [])
        ax.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
        ax.tick_params(axis="x", length=1.8, pad=.8, labelsize=8.5)
        if show_ticks:
            labels=ax.get_xticklabels()
            labels[0].set_horizontalalignment("left")
            labels[-1].set_horizontalalignment("right")
        for spine in ax.spines.values():
            spine.set_visible(False)


    for index, (treatment, title) in enumerate(panels):
        row, col = divmod(index, 2)
        x0, top = 3 + 89 * col, H - 7 - 94 * row
        text(x0, top, chr(65 + index), fontsize=10.5, weight="bold")
        text(x0 + 6, top, title, fontsize=9, weight="bold")
        group = curves[curves.exposure.eq(treatment)]
        wide = treatment == "exposure_crrt"
        # Use a shared scale within each treatment; never truncate a confidence band.
        limit = (.001, 2e6) if wide else (.30, 6.5)
        ticks = [.001, 1, 1e3, 1e6] if wide else [.5, 1, 2, 4]
        tick_labels = ["0.001", "1", "1,000", "1,000,000"] if wide else ["0.5", "1", "2", "4"]
        assert group.empty or (group.ci_low.min() > limit[0] and group.ci_high.max() < limit[1])
        for c, (database, db_label) in enumerate(cohorts):
            x, width = x0 + 12 + 39 * c, 33.5
            text(x, top - 8, db_label, weight="bold", fontsize=8.5)
            for r, (endpoint, time_label, metric) in enumerate(outcomes):
                y = top - 40 - 37 * r
                q = tests.loc[(treatment, database, endpoint)]
                this = group[group.dataset.eq(database) & group.outcome.eq(endpoint)]
                obs = support[support.exposure.eq(treatment) & support.dataset.eq(database)
                              & support.outcome.eq(endpoint)]
                ax = axes_at(x, y, width, 19)
                ax.set(xlim=(0, 1), ylim=limit, yscale="log", xticks=[])
                ax.set_yticks(ticks, labels=tick_labels if c == 0 else [])
                ax.yaxis.set_minor_locator(NullLocator())
                ax.tick_params(axis="y", labelsize=8.5, length=2, pad=.8)
                ax.grid(axis="y", lw=.5, color=light)
                ax.axhline(1, color=gray, ls=":", lw=.8)
                text(x, y + 28, f"{time_label} ({metric})", fontsize=8.5)
                if q.status == "estimated":
                    text(x, y + 22.8, f"Interaction P {p_format(q.interaction_p)}\nSlope-difference P {p_format(q.slope_difference_p)}",
                         color=gray, fontsize=7.5)
                    for k in (2, 3):
                        z = this[this.predicted_class.eq(k)].sort_values("probability")
                        if z.empty:
                            ax.text(.5, .22 if k == 2 else .08, f"Class {k}: no shared range",
                                    transform=ax.transAxes, ha="center", color=palette[k], fontsize=6.3)
                            continue
                        ax.fill_between(z.probability, z.ci_low, z.ci_high, color=palette[k], alpha=.15, lw=0)
                        ax.plot(z.probability, z.ratio, color=palette[k], ls=line_styles[k], lw=1.25)
                else:
                    assert this.empty
                    ax.text(.5, .56, "Not estimable", transform=ax.transAxes,
                            fontsize=7.2, weight="bold", ha="center")
                    if q.reason.startswith("Arm support"):
                        ax.text(.5, .26, "Sparse exposure arm", transform=ax.transAxes,
                                fontsize=6.6, ha="center", color=gray)
                support_strip(x, y - 4.7, width, obs, show_ticks=(r == 1))

    text(W / 2, 13, "Predicted class probability", fontsize=8.5, ha="center")
    fig.legend(handles=[Line2D([], [], color=palette[k], lw=1.45, ls=line_styles[k],
                               label=f"Class {k}") for k in (2, 3)],
               bbox_to_anchor=(.5, 6 / H), loc="center", ncols=2,
               frameon=False, fontsize=8.5, borderaxespad=0, handlelength=2)
    return fig
