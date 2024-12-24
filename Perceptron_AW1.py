import numpy as np
import matplotlib.pyplot as plt
from gc_contentNlength import fetch_data_for_multiple_queries #This is so I can access NCBI api for getting the desired information of genes, you may chose what organism and other query down below

class Perceptron():
    def __init__(self, input_size, learning_rate=0.8, epochs=1500): #I notied that to make any sort of boundary, we need to use large learning_rates and epochs, could be because the data is hardly linearly seperable and MAYBE something with how we norm the data
        self.weights = np.zeros(input_size)  # Initialize weights to zero
        self.bias = 0
        self.learning_rate = learning_rate
        self.epochs = epochs

    def act_func(self, x):  # Activation function: step function
        return 1 if x >= 0 else 0
    
    def prediction(self, input):  # Prediction: weighted sum + bias and activation
        output = np.dot(input, self.weights) + self.bias
        output = self.act_func(output)  # Apply activation function
        return output
    
    def normalize(self, X):
        # Normalize X by subtracting the mean and dividing by the standard deviation
        # This will scale the features to have zero mean and unit variance
        return (X - X.mean(axis=0)) / X.std(axis=0), X.mean(axis=0), X.std(axis=0)

    def fit(self, X, y):
        # Ensure X is a numpy array for proper slicing
        X = np.array(X)
        y = np.array(y)

        # Normalize the data before training and store the original mean and std
        X_normalized, X_mean, X_std = self.normalize(X)
        
        # Training process
        for epoch in range(self.epochs):
            for inputs, target in zip(X_normalized, y):
                prediction = self.prediction(inputs)
                error = target - prediction  # Calculate error
                # Update weights and bias
                for i in range(len(inputs)):
                    self.weights[i] += self.learning_rate * error * inputs[i]
                self.bias += self.learning_rate * error  # Update bias

        # After training, plot the decision boundary and data points
        self.plot_decision_boundary(X, y, X_mean, X_std)

    def plot_decision_boundary(self, X, y, X_mean, X_std):
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.set_xlim(X[:, 0].min() - 1, X[:, 0].max() + 1)
        ax.set_ylim(X[:, 1].min() - 1, X[:, 1].max() + 1)

        # Plot data points
        scatter = ax.scatter(X[:, 0], X[:, 1], c=y, cmap="coolwarm", edgecolors='k', marker='o', s=100)

        # Calculate decision boundary equation
        w1, w2 = self.weights
        b = self.bias

        # Inverse normalization to get boundary in original scale
        w1_original = w1 / X_std[0]
        w2_original = w2 / X_std[1]
        b_original = b - np.dot(w1_original * X_mean[0], w2_original * X_mean[1]) #Dont be confused that we tranform in back, we plot the original data points so we need a boundary that seperates them, not the transformed onces

        # Set a larger range for x_vals to simulate an "infinite" boundary
        x_vals = np.linspace(X[:, 0].min() - 5, X[:, 0].max() + 5, 100000)  # Larger range, you may add more zeros if you have a large and scatted inputs
        y_vals = -(w1_original * x_vals + b_original) / w2_original

        # Plot decision boundary (set color to red)
        ax.plot(x_vals, y_vals, 'r-', lw=2)

        # Display the decision boundary equation on the plot
        equation_text = f'Decision boundary: $y = -({w1_original:.2f}) x_1 - ({b_original:.2f})$ / {w2_original:.2f}'
        ax.text(0.05, 0.95, equation_text, transform=ax.transAxes, fontsize=12, verticalalignment='top', color='black')

        ax.set_title("Perceptron - Decision Boundary")
        ax.set_xlabel("GC Content")
        ax.set_ylabel("Sequence Length")

        plt.show()

# Parameters
size = 500
queries = [
    "Homo sapiens AND kinase",
    "hepatitis AND non-coding"
]
make_fasta = "N"
make_fasta_name = None

# Fetch data
X, y = fetch_data_for_multiple_queries(size, queries, make_fasta, make_fasta_name)

# Initialize and train perceptron
perceptron = Perceptron(input_size=2)
perceptron.fit(X, y)

# Test perceptron (optional, to see individual predictions)
for inputs in X:
    print(f"Input: {inputs}, Prediction: {perceptron.prediction(inputs)}")
