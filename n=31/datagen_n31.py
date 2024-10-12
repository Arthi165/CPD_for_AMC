# -*- coding: utf-8 -*-
"""

Generated in Colab.

"""

# The data generation of various AWGN

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Function to convert decimal number to binary vector
def decimal_to_binary_vector(decimal_no, binary_vector_length):
    if decimal_no < 0:
        raise ValueError("The decimal number must be non-negative.")
    binary_str = bin(decimal_no)[2:].zfill(binary_vector_length)
    if len(binary_str) > binary_vector_length:
        binary_str = binary_str[-binary_vector_length:]
    binary_vector = np.array([int(bit) for bit in binary_str])
    return binary_vector

# Function to calculate Hamming weight
def hamming_weight(vector):
    return np.sum(vector)

# Function to generate parity-check matrix from generator matrix
def generate_parity_check_matrix(G):
    k, n = G.shape  # Get dimensions of the generator matrix

    # Check if G is in standard form [I_k | P]
    I_k = np.eye(k, dtype=int)
    if not np.array_equal(G[:, :k], I_k):
        raise ValueError("Generator matrix G is not in standard form [I_k | P]")

    P = G[:, k:]  # Extract the P matrix
    P_T = P.T  # Transpose of P
    I_n_k = np.eye(n - k, dtype=int)  # Identity matrix of size (n - k)

    # Construct the parity-check matrix H
    H = np.hstack((P_T, I_n_k))
    return H

# Function to generate data for y vectors
def data_gen_channel_code_BPSK(G, SNR, N):
    k, n = G.shape  # Number of rows and columns in generator matrix

    y_vectors = []
    for _ in range(N):
        # Generate random message vector m of length k
        m = np.random.randint(2, size=k)

        # Encode message using generator matrix
        v = np.mod(np.dot(m, G), 2)

        # BPSK modulation
        v_prime = 1 - 2 * v

        # Calculate noise power sigma^2 from SNR
        snr_linear = 10 ** (SNR / 10.0)
        noise_power = 1 / snr_linear

        # Add AWGN noise
        noise = np.sqrt(noise_power / 2) * np.random.randn(*v_prime.shape)
        y = v_prime + noise

        y_vectors.append(y)

    y_vectors = np.vstack(y_vectors)
    return y_vectors

# Function to calculate probabilities, h_star vectors, and S_h set
def parity_check_probability(y_vectors, h_star, sigma_squared):
    # h_star as a 1 x n vector
    h_star = h_star.reshape(1, -1)

    # Set S_h: set of indices where h_star == 1
    S_h = {i for i, h_i in enumerate(h_star[0]) if h_i == 1}

    probabilities = []
    for y in y_vectors:
        tanh_values = np.tanh(y / sigma_squared)
        product_tanh = np.prod([tanh_values[j] for j in S_h])  # Product over indices in S_h
        probability = 0.5 + (0.5 * product_tanh)
        probabilities.append(probability)

    return np.array(probabilities)

def func1_h_in_C1_perp_notin_C2_perp(G1, G2):

    # Extract n1 and k1 from G1
    k1, n1 = G1.shape

    # Create h_star which is a vector containing all ones of order 1 x n1
    h_star = np.ones(n1, dtype=int)

    # Generate the parity-check matrix H1 from G1
    H1 = generate_parity_check_matrix(G1)

    valid_h_vectors = []
    min_weight = float('inf')

    # Iterate from 1 to 2^(n1 - k1)
    for i in range(1, 2**(n1 - k1)):
        m = decimal_to_binary_vector(i, n1 - k1)  # Convert decimal i to binary vector m
        h = np.mod(np.dot(m, H1), 2)  # Calculate the syndrome h = m * H1

        # Check the condition h * G2.T != 0
        if not np.all(np.mod(np.dot(h, G2.T), 2) == 0):
            valid_h_vectors.append(h)

            # Calculate the Hamming weight of h
            weight = hamming_weight(h)

            # Update h_star and min_weight if a lower weight is found
            if weight < min_weight:
                min_weight = weight
                h_star = h

    # Ensure the function returns valid outputs
    if h_star is None:
        raise ValueError("No valid syndrome vector found.")

    return h_star, valid_h_vectors, min_weight

try:
    # Define generator matrices G1 and G2
    import numpy as np

    # Define the generator matrix for (31, 21) code
    G2 = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0],
                   [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0],
                   [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0],
                   [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0],
                   [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0],
                   [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0],
                   [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

    G1 = np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1],
                   [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0],
                   [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0],
                   [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1],
                   [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
                   [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1],
                   [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
                   [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

    # Calculate h_star, valid_h_vectors, and min_weight for G1 and G2
    h_star_G1, valid_h_vectors_G1, min_weight_G1 = func1_h_in_C1_perp_notin_C2_perp(G1, G2)

    # Example y_vectors and sigma_squared
    N_train = 300000  # training size
    snr_value = 5  # SNR value
    
    # Generate data for G1
    y_vectors_G1 = data_gen_channel_code_BPSK(G1, snr_value, N_train)
    sigma_squared_G1 = 1 / (2 * (10 ** (snr_value / 10.0)))
    probabilities_G1 = parity_check_probability(y_vectors_G1, h_star_G1, sigma_squared_G1)

    # Save y_vectors_G1 to CSV
    df_y_vectors_G1 = pd.DataFrame(y_vectors_G1)
    df_y_vectors_G1.to_csv(f'csv_H1_no_awgn_SNR_{snr_value}dB.csv', index=False)
    print(f'y_vectors_G1 saved to csv_H1_no_awgn_SNR_{snr_value}dB.csv')

    # Generate data for G2
    y_vectors_G2 = data_gen_channel_code_BPSK(G2, snr_value, N_train)
    sigma_squared_G2 = 1 / (2 * (10 ** (snr_value / 10.0)))
    probabilities_G2 = parity_check_probability(y_vectors_G2, h_star_G1, sigma_squared_G2)

    # Save y_vectors_G2 to CSV
    df_y_vectors_G2 = pd.DataFrame(y_vectors_G2)
    df_y_vectors_G2.to_csv(f'csv_H2_no_awgn_SNR_{snr_value}dB.csv', index=False)
    print(f'y_vectors_G2 saved to csv_H2_no_awgn_SNR_{snr_value}dB.csv')

    # Generate parity-check matrices for G1 and G2
    H1 = generate_parity_check_matrix(G1)
    H2 = generate_parity_check_matrix(G2)

    # tanh is signalproc2
    # Save the results to CSV files
    df_G1 = pd.DataFrame(probabilities_G1, columns=['Probability'])
    df_G2 = pd.DataFrame(probabilities_G2, columns=['Probability'])

    df_G1.to_csv(f'csv_H1_tanh{snr_value}_SNR_{snr_value}dB.csv', index=False)
    df_G2.to_csv(f'csv_H2_tanh{snr_value}_SNR_{snr_value}dB.csv', index=False)

    print(f'Probabilities for G1 with SNR={snr_value}dB saved to csv_H1_tanh{snr_value}_SNR_{snr_value}dB.csv')
    print(f'Probabilities for G2 with SNR={snr_value}dB saved to csv_H2_tanh{snr_value}_SNR_{snr_value}dB.csv')

except ValueError as e:
    print(e)