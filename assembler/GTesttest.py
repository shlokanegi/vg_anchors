from assembler.gtest import GTest

# You can now create instances of your C++ class:
dummy_matrix = [[1, 20], [30, 1]]
gtest_instance = GTest(dummy_matrix, 0.01)

# And access its public members:
print(gtest_instance.success)
for hypothesis in gtest_instance.hypotheses:
    print(hypothesis.G, hypothesis.connectivityMatrix)
