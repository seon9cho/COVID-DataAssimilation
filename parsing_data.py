# ---------- Needed imports ----------  #
import pandas as pd

# ---------- Function for parsing infected data ----------  #
def parse_infected_data(file):
    """Parse weekly new case count data.
    
    All data is in 1-week increments here. 
    Testing data starts the week after the training data ends.

    Parameters:
        file (str): The path to "Summary_Case Counts Per Week.csv"

    Returns: 
        train (np.ndarray) - An array of weekly new COVID case counts 
            for training (the weeks from 3/8/2020 to 12/12/2021).
            Use for your gradient descent algorithm and modeling.
        test (np.ndarray) - An array of weekly new COVID case counts
            for testing (the weeks from 12/19/2021 to 1/19/2022). 
            Only use to test the performance of your model on its 
            learned parameters after you run your gradient descent 
            algorithm.
        t_train (list) - A list of numbered weeks
            corresponding to the training data (starting with 1) 
        t_test (list) - A list of numbered weeks 
            corresponding to the testing data (starting 1 after
            the training data)
    """
    # Get the case count data from the csv into Python
    data = pd.read_csv(file, index_col=0)

    # Split the data into a training set(before 12/12/2021) 
    # and a testing set (12/19/2021-1/19/2022)
    train = data[:'2021-12-12']['New Cases'].to_numpy()
    test = data['2021-12-19':'2022-01-19']['New Cases'].to_numpy()

    # Get week enumerations for this data (for ease of use) 
    t_train = [i for i in range(1, len(train) + 1)]
    t_test = [i for i in range(len(train)+1, len(test)+len(train)+1)]
    return train, test, t_train, t_test



# ---------- Example for how to use the above function ----------  #
csv_file = 'COVID_Data_SLC/Summary_Case Counts per Week.csv'
train, test, t_train, t_test = parse_infected_data(csv_file)



