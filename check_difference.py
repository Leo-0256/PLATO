def check_difference(filepath1, filepath2):
    """
    Check the [cost] difference of two path file.
    """
    f1 = open(filepath1, "r")
    f2 = open(filepath2, "r")
    cost1 = []
    cost2 = []
    for row in f1:
        cost1.append(int(row.strip().split()[-1]))
    for row in f2:
        cost2.append(int(row.strip().split()[-1]))
    if len(cost1) != len(cost2):
        print("FALSE")
        return
    for i in range(len(cost1)):
        if cost1[i] != cost2[i]:
            print("FALSE")
            return
    print("TRUE")


if __name__ == "__main__":
    check_difference("./PATH/demo_2x2_path.txt", "./PATH/demo_2x2_nx_path.txt")
    check_difference("./PATH/demo_3x3_path.txt", "./PATH/demo_3x3_nx_path.txt")
    check_difference("./PATH/demo_4x4_path.txt", "./PATH/demo_4x4_nx_path.txt")
    check_difference("./PATH/demo_5x5_path.txt", "./PATH/demo_5x5_nx_path.txt")
    check_difference("./PATH/demo_10x10_path.txt", "./PATH/demo_10x10_nx_path.txt")
    check_difference("./PATH/demo_16x16_path.txt", "./PATH/demo_16x16_nx_path.txt")
    check_difference("./PATH/demo_25x25_path.txt", "./PATH/demo_25x25_nx_path.txt")
