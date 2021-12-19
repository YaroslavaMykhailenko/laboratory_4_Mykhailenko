import numpy as np
from matplotlib import pyplot as plt

def AbsorptionMarkov():

    # 1.1 Transition matrix

    P_matrix = [[1, 0, 0, 0], [0.22, 0.11, 0.58, 0.09], [0.15, 0.27, 0.20, 0.38], [0.19, 0.57, 0.19, 0.05]]
    vector = [0.1, 0.3, 0.2, 0.4]

        

    def printMatrix ( matrix ): 
        for i in range ( len(matrix) ): 
            for j in range ( len(matrix[i]) ): 
                print(matrix[i][j], end="\t")
            print ()

    print("Transition matrix:")
    printMatrix(P_matrix)


    # Model a Markov chain with absorption. Quantity implementations - more than 100

    number_of_realization = 100
    arr_realization = []
    n = 30

    while len(arr_realization) < number_of_realization:
        i = vector.index(np.random.choice(vector, 1, p = vector)[0])
        realization = []

        while i != 0:
            realization.append(i)
            vector = P_matrix[i]
            i = vector.index(np.random.choice(vector, 1, p = vector)[0])
        realization.append(i)

        if len(realization) > 0:
            arr_realization.append(realization)


    plt.figure(figsize=(10, 5))
    for i in range(n):
        plt.plot([ j for j in range(len(arr_realization[i])) ], arr_realization[i])
    plt.show()


    # Estimates of the transition matrix

    P_matrix_realization = np.zeros_like(P_matrix)
    for realization in arr_realization:
        matrix_tmp = np.zeros_like(P_matrix)
        matrix_tmp[0][0] = 1

        for (i, j) in zip(realization, realization[1:]):
            matrix_tmp[i][j] += 1
        P_matrix_realization = P_matrix_realization + matrix_tmp

    for index, row in enumerate(P_matrix_realization):
        P_matrix_realization[index] = row / sum(row)


    print("P_matrix of received realizations : \n",P_matrix_realization)


     # 1.2 Fundamental_matrix

    Q = []

    for i in range(1, len(P_matrix_realization)):
        arr_row = []

        for j in range(1,len(P_matrix_realization)):
            arr_row.append(P_matrix_realization[i][j])

        Q.append(arr_row)


    fundamental_matrix = np.linalg.inv(np.identity(len(Q)) - Q)
    print("Fundamental_matrix: \n", fundamental_matrix)



    # 1.3 The average number of steps that the chain is in the state j when the process started with the state 

    condition_n = int(input("Start condition: "))
    condition_j = int(input("Last condition: "))
    print("The average number of steps that the chain is in the state j when the process started with the state n: \n", fundamental_matrix[condition_n - 2][condition_j - 2], "\n")


    # 1.4 Average absorption time

    avarage_absorption = np.dot(fundamental_matrix, np.ones((len(fundamental_matrix))))
    print("Аverage step numbers:", avarage_absorption)



    # 1.5 Probability of absorption

    probability = []

    for i in range(1, len(P_matrix_realization)):
        row = []

        for j in range(len(P_matrix_realization) - 3):
            row.append(P_matrix_realization[i][j])

        probability.append(row)


    B = np.dot(fundamental_matrix, probability)

    print("Probability of absorption:\n", B)



def RegularMarkov():

    # 1.1 Transition matrix

    P_matrix = [[0.36, 0.26, 0.38], [0.93, 0.01, 0.06], [0.05, 0.54, 0.41]]  
    vector = [0.25, 0.26, 0.49]

    def printMatrix ( matrix ): 
        for i in range ( len(matrix) ): 
            for j in range ( len(matrix[i]) ): 
                print(matrix[i][j], end="\t")
            print ()

    print("Transition matrix:")
    printMatrix(P_matrix)


    # Model a regular Markov chain. Quantity implementations - more than 100

    number_of_realization = 100
    arr_realization = []
    n = 30

    while len(arr_realization) < number_of_realization:
        i = vector.index(np.random.choice(vector, 1, p=vector)[0])
        realization = []
        while len(realization) < 10:
            realization.append(i)
            vector = P_matrix[i]
            i = vector.index(np.random.choice(vector, 1, p=vector)[0])
        realization.append(i)
        arr_realization.append(realization)
    

    plt.figure(figsize=(15, 5))

    for i in range(n):
        plt.plot([j for j in range(len(arr_realization[i]))], arr_realization[i])

    plt.show()


    # Estimates of the transition matrix

    P_matrix_realization = np.zeros_like(P_matrix)
    for realization in arr_realization:
        matrix_tmp = np.zeros_like(P_matrix)

        for (i, j) in zip(realization, realization[1:]):
            matrix_tmp[i][j] += 1
        P_matrix_realization = P_matrix_realization + matrix_tmp

    for index, row in enumerate(P_matrix_realization):
        P_matrix_realization[index] = row / sum(row)

    print("P_matrix of received realizations : \n",P_matrix_realization)



    # 2.1 Final probabilities

    w = vector
    for i in range(1000):
        w = np.dot(w, P_matrix_realization).tolist()

    probability = []

    for J in range(len(P_matrix)):
        probability.append(w)

    probability = np.array(probability)

    print('Final matrix:')
    print(probability)


    # 1.2 Fundamental_matrix

    fundamental_matrix = np.linalg.inv(np.identity(len(P_matrix)) - P_matrix_realization + probability)
    print("Fundamental_matrix: \n", fundamental_matrix)


    #  2.3 average time in a given state for n = 4 steps

    fundamential_dg = np.diagflat(np.diag(fundamental_matrix))

    l = len(fundamental_matrix)
    E = [[1] * (l) for i in range(l)]

    M = np.dot((np.eye(l) - fundamental_matrix + np.dot(E, fundamential_dg)), np.diagflat(1 / probability[0]))

    print('Fverage time in a given state for n = 4 steps:')
    print(np.dot(vector, fundamental_matrix) - probability[0] + 4 * probability[0], "\n")


    # 2.4 the average time of the circuit in a given state (in the state j when the process started with the state n);

    condition_n = int(input("Start condition: "))
    condition_j = int(input("Last condition: "))
    print("The average time of the circuit in a given state (in the state j when the process started with the state n) ",M[ condition_n-1][condition_j-1], "\n")


    # 2.5 the average time of exit of the circuit to a given state in stationary mode (when the initial state is not set).

    print("The average time of exit of the circuit to a given state in stationary mode (when the initial state is not set): ",np.dot(probability[0], M))







def start():
    print("Choose task you want to solve: \n1 - Absorption Markov Chain\n2 - Regular Markov Chain")
    number = int(input(""))

    if number == 1:
        AbsorptionMarkov()
        again()
    elif number == 2:
        RegularMarkov()
        again()   
    else:
        print("Wrong number!")


def again():
    print("Would you like to try again?:\n1 - yes\n2 - no")
    choice = int(input(""))
    if choice == 1:
        start()
    
    elif choice == 2:
        return 0 


start()