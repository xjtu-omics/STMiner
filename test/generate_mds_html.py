import plotly
# for data visualization
import plotly.express as px
# for creating a swiss roll
from sklearn.datasets import make_swiss_roll

# Make a swiss roll
X, y = make_swiss_roll(n_samples=2000, noise=0.05)
# Make it thinner
X[:, 1] *= .5

# Create a 3D scatter plot
fig = px.scatter_3d(None, x=X[:, 0], y=X[:, 1], z=X[:, 2], color=y, )

# Update chart looks
fig.update_layout(
    showlegend=False,
    scene_camera=dict(up=dict(x=0, y=0, z=1),
                      center=dict(x=0, y=0, z=-0.1),
                      eye=dict(x=1.25, y=1.5, z=1)),
    margin=dict(l=0, r=0, b=0, t=0),
    scene=dict(xaxis=dict(backgroundcolor='white',
                          color='black',
                          gridcolor='#f0f0f0',
                          title_font=dict(size=10),
                          tickfont=dict(size=10),
                          ),
               yaxis=dict(backgroundcolor='white',
                          color='black',
                          gridcolor='#f0f0f0',
                          title_font=dict(size=10),
                          tickfont=dict(size=10),
                          ),
               zaxis=dict(backgroundcolor='lightgrey',
                          color='black',
                          gridcolor='#f0f0f0',
                          title_font=dict(size=10),
                          tickfont=dict(size=10),
                          )))

# Update marker size
fig.update_traces(marker=dict(size=3, line=dict(color='black', width=0.1)))

fig.update(layout_coloraxis_showscale=False)
fig.update_layout(width=400, height=300)

plotly.offline.plot(fig, filename='mds.html', image_width=400, image_height=300)
