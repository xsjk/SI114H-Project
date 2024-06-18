import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt

# Define the neural network
class PINN(nn.Module):
    def __init__(self, layers_v, layers_q):
        super(PINN, self).__init__()
        self.activation = nn.Tanh()
        
        self.layers_v = nn.ModuleList()
        for i in range(len(layers_v) - 1):
            self.layers_v.append(nn.Linear(layers_v[i], layers_v[i + 1]))

        self.layers_q = nn.ModuleList()
        for i in range(len(layers_q) - 1):
            self.layers_q.append(nn.Linear(layers_q[i], layers_q[i + 1]))

    def forward(self, x, t):
        inputs = torch.cat([x, t], dim=1)
        
        v = inputs
        for i in range(len(self.layers_v) - 1):
            v = self.activation(self.layers_v[i](v))
        v = self.layers_v[-1](v)
        
        q = inputs
        for i in range(len(self.layers_q) - 1):
            q = self.activation(self.layers_q[i](q))
        q = self.layers_q[-1](q)
        
        return v, q

# Define the load function f(x, t) in the differential equation
def f(x, t):
    return torch.sin(np.pi * x) * torch.cos(np.pi * t)

# Define the loss function
def loss_function(model, x, t, f):
    v, q = model(x, t)
    v_t = torch.autograd.grad(v, t, grad_outputs=torch.ones_like(v), create_graph=True)[0]
    v_x = torch.autograd.grad(v, x, grad_outputs=torch.ones_like(v), create_graph=True)[0]
    q_x = torch.autograd.grad(q, x, grad_outputs=torch.ones_like(q), create_graph=True)[0]
    q_t = torch.autograd.grad(q, t, grad_outputs=torch.ones_like(q), create_graph=True)[0]
    
    loss_region1 = torch.mean((0.00525 * v_t - q_x - f(x, t))**2 * (x > 0.65).float())
    loss_region2 = torch.mean((q_t - 60 * v_x)**2 * (x <= 0.65).float())
    
    loss = loss_region1 + loss_region2
    return loss

# Training the model
def train(model, x, t, f, epochs, lr):
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    for epoch in range(epochs):
        optimizer.zero_grad()
        loss = loss_function(model, x, t, f)
        loss.backward()
        optimizer.step()
        if epoch % 100 == 0:
            print(f'Epoch {epoch}, Loss: {loss.item()}')
    return model

# Define the neural network structure
layers_v = [2, 20, 20, 20, 1]
layers_q = [2, 20, 20, 20, 1]
model = PINN(layers_v, layers_q)

# Generate training data
x = torch.linspace(0, 1, 100).view(-1, 1)
t = torch.linspace(0, 1, 100).view(-1, 1)
x, t = torch.meshgrid(x.squeeze(), t.squeeze())
x = x.reshape(-1, 1)
t = t.reshape(-1, 1)
x.requires_grad = True
t.requires_grad = True

# Train the model
epochs = 1000
lr = 0.001
model = train(model, x, t, f, epochs, lr)

# Plot the results for v
x_test = torch.linspace(0, 1, 100).view(-1, 1)
t_test = torch.linspace(0, 1, 100).view(-1, 1)
x_test, t_test = torch.meshgrid(x_test.squeeze(), t_test.squeeze())
x_test = x_test.reshape(-1, 1)
t_test = t_test.reshape(-1, 1)
v_pred, q_pred = model(x_test, t_test)

v_pred = v_pred.detach().numpy().reshape(100, 100)
q_pred = q_pred.detach().numpy().reshape(100, 100)
x_test = x_test.detach().numpy().reshape(100, 100)
t_test = t_test.detach().numpy().reshape(100, 100)

#plot v and q

plt.figure()
plt.subplot(1, 2, 1)
plt.pcolor(x_test, t_test, v_pred, cmap='jet')
plt.colorbar()
plt.xlabel('x')
plt.ylabel('t')
plt.title('v')
plt.subplot(1, 2, 2)
plt.pcolor(x_test, t_test, q_pred, cmap='jet')
plt.colorbar()
plt.xlabel('x')
plt.ylabel('t')
plt.title('q')

plt.show()
