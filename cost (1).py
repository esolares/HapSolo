import start
def cost(find_count_dict):
    missingB = find_count_dict["missing"]
    dupeB = find_count_dict["duplicate"]
    fragB = find_count_dict["fragmented"]
    compB = find_count_dict["complete"]
    cost = (start.weight_missing*missingB)+(start.weight_duplicate*dupeB)+(fragB*start.weight_fragmented)
    if (start.weight_single*compB) != 0:
        cost = cost/(start.weight_single*compB)
    return cost