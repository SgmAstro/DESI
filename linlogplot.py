import matplotlib.pyplot as plt

# TODO: change plot to scatter.

def linlogplot(x,y):
    fig, axs = plt.subplots(2, 2, gridspec_kw = {'wspace':0, 'hspace':0}, figsize=(10, 10))
    axs[0, 0].scatter(x, y, s=0.25)
    axs[0, 0].set_yscale('log')
    axs[0, 0].set_xlim(0, 1)
    axs[0, 0].set_ylim(1.00001, 10)
    axs[0, 0].xaxis.set_visible(False)

    axs[0, 1].scatter(x, y, 'tab:orange')
    axs[0, 1].set_xlim(1.00001, 10)
    axs[0, 1].set_ylim(1.00001, 10)
    axs[0, 1].set_xscale('log')
    axs[0, 1].set_yscale('log')
    axs[0, 1].spines['left'].set_visible(False)
    axs[0, 1].xaxis.set_visible(False)
    axs[0, 1].yaxis.set_visible(False)

    axs[1, 0].scatter(x, y, 'tab:green')
    axs[1, 0].set_xlim(0, 1)
    axs[1, 0].set_ylim(0, 1)
    axs[1, 0].spines['top'].set_visible(False)

    axs[1, 1].scatter(x, y, 'tab:red')
    axs[1, 1].set_xscale('log')
    axs[1, 1].set_xlim(1.00001, 10)
    axs[1, 1].set_ylim(0, 1)
    axs[1, 1].spines['left'].set_visible(False)
    axs[1, 1].spines['top'].set_visible(False)
    axs[1, 1].set_yticks([])