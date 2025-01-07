import torch
import matplotlib.pyplot as plt

# Grid parameters
num_rows, num_cols = 10, 10
branch_prob = 1

# State of the grid: single active neuron in the middle of column 0
states = torch.zeros((num_rows, num_cols))
states[num_rows // 2, 0] = 1

# Feed-forward activation to the next column only
for c in range(1, num_cols):
    for r in range(num_rows):
        if states[r, c - 1] == 1:
            if torch.rand(1) < branch_prob:
                states[r, c] = 1

# Plot the final pattern
plt.imshow(states, cmap='Greys', origin='lower')
plt.title('Single-Layer Feed-Forward Activation')
plt.colorbar(label='Active (1) or Inactive (0)')
plt.show()
