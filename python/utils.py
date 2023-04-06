import os
import pickle

def save_results_to_date_file(data_to_save, file_name: str = None):
    from datetime import datetime, date

    today = date.today()
    d1 = today.strftime("%d_%m_%Y")

    if not file_name:
        now = datetime.now()
        file_name = now.strftime("%H:%M")

    path = os.path.join('angleAction','saved_results', d1)
    if not os.path.isdir(path):
        os.makedirs(path)

    path = os.path.join(path, file_name)
    with open(path + '.pickle', 'wb') as handle:
        pickle.dump(data_to_save, handle, protocol=pickle.HIGHEST_PROTOCOL)

