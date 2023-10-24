def merge_isolands(isolands):
    """
    Merge isolands

    Parameters
    --------------
    isolands : a list of number contain isolands need to merge
        A list in form like [[1, 3], [2, 4], [5, 9]]

    Returns
    ----------
    isolands_merged : a list of number contain merged isolands
    """
    
    labeled_border = []

    index = 0
    for isol in isolands:
        labeled_border.append((isol[0], "l", index))
        labeled_border.append((isol[1], "r", index))
        index += 1

    labeled_border += 1

    left = 0
    right = 0
    isolands_merged = []
    isolands_merged.append(labeled_border[0][0])
    one_isol = []
    for border in labeled_border:
        if border[1] == "l":
            left += 1
        elif border[1] == "r":
            right += 1
        if left - right == 1:
            one_isol.append(border[0])
        if left == right:
            one_isol.append(border[0])
            isolands_merged.append(one_isol)
            one_isol = []

    return isolands_merged


def change_coordinate(left, positions, cor="relative"):
    new_pos = []
    if cor == "relative":
        for pos in positions:
            new_pos.append([pos[0] - left + 1, pos[1] - left + 1])
    elif cor == "absolute":
        for pos in positions:
            new_pos.append([pos[0] + left - 1, pos[1] + left - 1])

    return new_pos
