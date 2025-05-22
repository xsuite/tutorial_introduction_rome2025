# import useful libraries
import numpy as np  # arrays and math
import matplotlib.pyplot as plt  # plots
import xtrack as xt  # tracking module of Xsuite
import matplotlib.gridspec as gridspec

from matplotlib.widgets import Slider


def track_and_plot_eigen(line, part, fignumber=1):

    line.track(
        part.copy(), turn_by_turn_monitor="ONE_TURN_EBE"
    )  # collect element by element data
    data = line.record_last_track
    m = FODOline.twiss4d(betx=1, bety=1).get_R_matrix(
        start="start", end="stop"
    )  # transfer matrix
    eig = np.linalg.eigvals(m[:4, :4])

    # Checking the compatibility of the line with prebuilt kernels
    # line.tracker.check_compatibility_with_prebuilt_kernels()  # check the compatibility of the line

    LongFODOline = 10 * FODOline
    LongFODOline.track(
        part.copy(), turn_by_turn_monitor="ONE_TURN_EBE"
    )  # collect element by element data
    dataLONG = LongFODOline.record_last_track

    # Create the overall figure and outer GridSpec (3 rows total)
    fig = plt.figure(figsize=(8, 8), num=fignumber)
    gs_main = gridspec.GridSpec(3, 1, height_ratios=[2, 1, 1], figure=fig)

    # Top block: 2x2 with spanning on right
    gs_top = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs_main[0])

    HorTrackingCELL = fig.add_subplot(gs_top[0, 0])  # Top-left
    VerTrackingCELL = fig.add_subplot(gs_top[1, 0])  # Bottom-left
    EigenValueCELL = fig.add_subplot(gs_top[:, 1])  # Right column spans both rows

    # Bottom rows
    HorTrackingLINE = fig.add_subplot(gs_main[1])  # Full-width row
    VerTrackingLINE = fig.add_subplot(gs_main[2])  # Full-width row

    # Label axes for visualization
    HorTrackingCELL.set_title("FODO cell tracking")
    # VerTrackingCELL.set_title("FODO cell: y track")
    EigenValueCELL.set_title("FODO cell R-matrix eigenvalues")
    HorTrackingLINE.set_title("10 FODO cells tracking")
    # VerTrackingLINE.set_title("10 FODO cells: y track")

    # Plotting the FODO cell tracking
    HorTrackingCELL.plot(data.s.T, data.x.T, "-")
    VerTrackingCELL.plot(data.s.T, data.y.T, "-")
    VerTrackingCELL.set_xlabel("s [m]")
    HorTrackingCELL.set_ylabel("x [m]")
    VerTrackingCELL.set_ylabel("y [m]")
    HorTrackingCELL.fill_between(
        [9.5, 10.2], [-0.08, -0.08], [0.08, 0.08], color="b", alpha=0.5
    )
    VerTrackingCELL.fill_between(
        [9.5, 10.2], [-0.08, -0.08], [0.08, 0.08], color="b", alpha=0.5
    )

    # Plotting the eigenvalues
    (EigenValuePlot1,) = EigenValueCELL.plot(eig[:2].real, eig[:2].imag, "ok")
    (EigenValuePlot2,) = EigenValueCELL.plot(eig[2:].real, eig[2:].imag, "or")
    EigenValueCELL.grid()
    EigenValueCELL.set_xlim(-6, 6)
    EigenValueCELL.set_ylim(-1.5, 1.5)
    # Plotting the 10 FODO cells tracking
    HorTrackingLINE.plot(dataLONG.s.T, dataLONG.x.T, "-")
    VerTrackingLINE.plot(dataLONG.s.T, dataLONG.y.T, "-")
    VerTrackingLINE.set_xlabel("s [m]")
    HorTrackingLINE.set_ylabel("x [m]")
    VerTrackingLINE.set_ylabel("y [m]")

    # Set the y-limits for the tracking plots
    [ax.set_ylim(-0.2, 0.2) for ax in [HorTrackingLINE, VerTrackingLINE]]
    # Set the y-limits for the cell tracking plots
    [ax.set_ylim(-0.2, 0.2) for ax in [HorTrackingCELL, VerTrackingCELL]]

    plt.tight_layout()

    return (
        fig,
        HorTrackingCELL,
        VerTrackingCELL,
        EigenValuePlot1,
        EigenValuePlot2,
        HorTrackingLINE,
        VerTrackingLINE,
    )


def updateFull(val):
    kf = slider_kf_Full.val  # Get the current beta value
    kd = slider_kd_Full.val  # Get the current theta value

    env["kf"] = kf  # Update the beta value in the environment
    env["kd"] = kd  # Update the theta value in the environment
    # print(f"kf={kf} 1/m, kd={kd} 1/m") # for debugging

    try:
        m = FODOline.twiss4d(betx=1, bety=1).get_R_matrix(
            start="start", end="stop"
        )  # transfer matrix
        eig = np.linalg.eigvals(m[:4, :4])

        FODOline.track(
            p1.copy(), turn_by_turn_monitor="ONE_TURN_EBE"
        )  # collect element by element data
        data = FODOline.record_last_track

        LongFODOline = 10 * FODOline
        LongFODOline.track(
            p1.copy(), turn_by_turn_monitor="ONE_TURN_EBE"
        )  # collect element by element data
        dataLONG = LongFODOline.record_last_track

        # print(data.x[0]) # Fog debugging

        lines = HorTrackCELL.get_lines()
        for ii in range(len(lines)):
            line = lines[ii]
            line.set_ydata(data.x[ii])  # Update the horizontal tracking plot
            # Get the line object

        lines = VerTrackCELL.get_lines()
        for ii in range(len(lines)):
            line = lines[ii]
            line.set_ydata(data.y[ii])  # Update the vertical tracking plot
            # Get the line object

        lines = HorTrackLINE.get_lines()
        for ii in range(len(lines)):
            line = lines[ii]
            line.set_ydata(dataLONG.x[ii])  # Update the horizontal tracking plot
            # Get the line object

        lines = VerTrackLINE.get_lines()
        for ii in range(len(lines)):
            line = lines[ii]
            line.set_ydata(dataLONG.y[ii])  # Update the horizontal tracking plot
            # Get the line object

        EigenCELLplot1.set_data(eig[:2].real, eig[:2].imag)  # Update the  plot
        EigenCELLplot2.set_data(eig[2:].real, eig[2:].imag)  # Update the  plot

    except Exception as e:
        print(f"Error during tracking: {e}")

    WholeFig.canvas.draw_idle()  # Redraw the figure


## Build the trasport line (FODO line)
# define a FODO line
env = xt.Environment()  # create a simulation environment
env.new("start", xt.Marker)  # create a marker
env.new("stop", xt.Marker)  # create a marker
FODOline = env.new_line(
    name="FODOline", components=[env.place("start", at=0), env.place("stop", at=10)]
)

# FODOline.config.XTRACK_USE_EXACT_DRIFTS = True  # use exact drifts

# Define and place quadrupoles
# Starting values for quadrupole strengths
env["kf"] = 2.8
env["kd"] = -2.35
print(f"Starting values: kf={env['kf']} 1/m, kd={env['kd']} 1/m")

if "qd" not in env.element_dict:
    env.new(
        "qd", xt.Quadrupole, length=0.20, k1="kd"
    )  # k1 is scaled gradient k1=e/q G [1/m]
if "qf" not in env.element_dict:
    env.new(
        "qf", xt.Quadrupole, length=0.20, k1="kf"
    )  # k1 is scaled gradient k1=e/p G [1/m]

FODOline.insert("qf", at=3)
FODOline.insert("qd", at=6)

# Check the line: table and survey plot
# Check the line
# FODOline.get_table()
# FODOline.get_table().cols['name s element_type isthick'] # get the table of elements
# FODOline.survey().plot()  # plot the line
# plt.show()

# Define the particles and track them
# define particles
part = xt.Particles(mass0=xt.PROTON_MASS_EV, p0c=1e9, x=np.linspace(-0.001, 0.001, 11))
p1 = part.copy()  # copy the particles
p1.x = np.random.uniform(-0.01, 0.01, 11)  # random the horizontal position
p1.y = np.random.uniform(-0.01, 0.01, 11)  # random the vertical position
p1.px = np.random.uniform(-0.01, 0.01, 11)  # random the horizontal momentum
p1.py = np.random.uniform(-0.01, 0.01, 11)  # random the vertical momentum

FODOline.particle_ref = xt.Particles(mass0=xt.PROTON_MASS_EV, p0c=1e9)


# Define the figure and axis
(
    WholeFig,
    HorTrackCELL,
    VerTrackCELL,
    EigenCELLplot1,
    EigenCELLplot2,
    HorTrackLINE,
    VerTrackLINE,
) = track_and_plot_eigen(FODOline, p1)
WholeFig.subplots_adjust(bottom=0.3)

# Create the slider axes and sliders
ax_slider_kf_Full = WholeFig.add_axes(
    [0.2, 0.17, 0.65, 0.03]
)  # [left, bottom, width, height]
ax_slider_kd_Full = WholeFig.add_axes(
    [0.2, 0.14, 0.65, 0.03]
)  # [left, bottom, width, height]

slider_kf_Full = Slider(
    ax_slider_kf_Full, "k focusing", 0, 5.0, valinit=env["kf"], valstep=0.01
)
slider_kd_Full = Slider(
    ax_slider_kd_Full, "k defocusion", -5, 0, valinit=env["kd"], valstep=0.01
)

slider_kf_Full.on_changed(updateFull)
slider_kd_Full.on_changed(updateFull)

plt.show()
