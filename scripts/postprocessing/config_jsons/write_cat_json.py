import argparse
import json

def generate_categories(pt_boundaries, isData=False):
    if isData: sigma_m_over_m_string = "sigma_m_over_m_smeared_decorr"
    else: sigma_m_over_m_string = "sigma_m_over_m_corr_smeared_decorr"
    categories = {}
    mass_resolution_categories = ['cat0', 'cat1', 'cat2']
    mass_resolution_thresholds = [0.010, 0.014]
    
    for i in range(len(pt_boundaries) - 1):
        for j, mass_resolution_cat in enumerate(mass_resolution_categories):
            category_name = f"PTH_{str(pt_boundaries[i]).replace('.', 'p')}_{str(pt_boundaries[i + 1]).replace('.', 'p')}_{mass_resolution_cat}"
            category_filters = [
                ["lead_corr_mvaID_run3", ">", 0.25],
                ["sublead_corr_mvaID_run3", ">", 0.25],
                ["pt", ">", pt_boundaries[i]],
                ["pt", "<", pt_boundaries[i + 1]]
            ]
            if j == 0:
                category_filters.append([sigma_m_over_m_string, "<", mass_resolution_thresholds[j]])
            elif j == len(mass_resolution_categories) - 1:
                category_filters.append([sigma_m_over_m_string, ">", mass_resolution_thresholds[j - 1]])
            else:
                category_filters.append([sigma_m_over_m_string, ">", mass_resolution_thresholds[j - 1]])
                category_filters.append([sigma_m_over_m_string, "<", mass_resolution_thresholds[j]])
            
            categories[category_name] = {
                "cat_filter": category_filters
            }
    return categories

def save_to_json(output_file, categories):
    data = {"cat_dict": categories}
    with open(output_file, 'w') as outfile:
        json.dump(data, outfile, indent=4)

def main():
    parser = argparse.ArgumentParser(description='Generate JSON file with specified PT boundaries')
    parser.add_argument('output_file', help='Output file name and location')
    parser.add_argument('pt_boundaries', nargs='+', type=float, help='PT boundaries list')
    parser.add_argument('--isData', action="store_true", default=False, help="Add this flag if you are running over data, this changes the name of the sigma_m/m variable that is read.")

    args = parser.parse_args()

    print(args.isData)

    categories = generate_categories(args.pt_boundaries, args.isData)
    save_to_json(args.output_file, categories)

if __name__ == "__main__":
    main()