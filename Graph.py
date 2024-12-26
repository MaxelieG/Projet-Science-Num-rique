import matplotlib.pyplot as plt

def heatmap_plot (figure, temperature_matrix, title, position):

    axes = figure.add_subplot(1, 3, position)
    axes.matshow(temperature_matrix, cmap = "hot")
    axes.set_title(title)
    axes.set_xlabel('Lx')
    if position == 1:
        axes.set_ylabel('Ly')
    figure.colorbar(axes.matshow(temperature_matrix, cmap = "hot"), ax=axes, orientation = "horizontal", fraction=0.046, pad=0.04)
    
    return None