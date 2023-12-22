import os
import pathlib

class cleanup():
    '''
    Folder(n_layers) we would want to clean: 1_full_sequence_dataset(1), 2_Cut_Data(2), 3_cut_dataset(3),4_relabeled_dataset(3),
    5_inferred(2)
    Be cautious with 6_results(1), 7_evals(2), logs(2)
    '''
    def clean(self, folder, n_layers):
        if n_layers == 1:
            for file in os.listdir(folder):
                os.remove(os.path.join(folder, file))
            os.rmdir(folder)
        elif n_layers == 2:
            for folder2 in os.listdir(folder):
                for file in os.listdir(os.path.join(folder, folder2)):
                    os.remove(os.path.join(folder, folder2, file))
                os.rmdir(os.path.join(folder, folder2))
            os.rmdir(folder)
        elif n_layers == 3:
            for folder2 in os.listdir(folder):
                for folder3 in os.listdir(os.path.join(folder, folder2)):
                    for file in os.listdir(os.path.join(folder, folder2, folder3)):
                       os.remove(os.path.join(folder, folder2, folder3, file))
                    os.rmdir(os.path.join(folder, folder2, folder3))
                os.rmdir(os.path.join(folder, folder2))
            os.rmdir(folder)
    '''
    removes all logs
    '''
    def clean_logs(self):
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        for folder1 in os.listdir(os.path.join(py_eval_path, 'logs')):
            for folder2 in os.listdir(os.path.join(py_eval_path, 'logs', folder1)):
                for folder3 in os.listdir(os.path.join(py_eval_path, 'logs', folder1, folder2)):
                    for file in os.listdir(os.path.join(py_eval_path, 'logs', folder1, folder2, folder3)):
                        os.remove(os.path.join(py_eval_path, 'logs', folder1, folder2, folder3, file))
                    os.rmdir(os.path.join(py_eval_path, 'logs', folder1, folder2, folder3))
                os.rmdir(os.path.join(py_eval_path, 'logs', folder1, folder2))
            os.rmdir(os.path.join(py_eval_path, 'logs', folder1))
        os.rmdir(os.path.join(py_eval_path, 'logs'))
    '''
    checks if folder needs to be deleted
    '''
    def check_existence(self, folders):
        exist = True
        for folder in folders:
            if not os.path.exists(folder):
                exist = False
                break
        return exist

    '''
    Deletes all files generated during training
    '''
    def post_training_cleanup(self):
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        check = self.check_existence([os.path.join(py_eval_path, '3_cut_dataset'), os.path.join(py_eval_path, '4_relabeled_dataset')])
        if check:
            self.clean(os.path.join(py_eval_path, '3_cut_dataset'), 3)
            self.clean(os.path.join(py_eval_path, '4_relabeled_dataset'), 3)
        check2 = self.check_existence([os.path.join(py_eval_path, '1_full_sequence_dataset'), os.path.join(py_eval_path, '0_data_mat')])
        if check2:
            self.clean(os.path.join(py_eval_path, '1_full_sequence_dataset'), 1)
            self.clean(os.path.join(py_eval_path, '0_data_mat'), 1)
        check3 = self.check_existence([os.path.join(py_eval_path, 'logs')])
        if check3:
            self.clean_logs()

    '''
    Deletes all files generated during evaluation
    '''
    def post_inference_cleanup(self):
        py_eval_path = str(pathlib.Path(__file__).parent.parent)
        paths = [os.path.join(py_eval_path, '2_Cut_Data'), os.path.join(py_eval_path, '5_inferred'), os.path.join(py_eval_path, '6_results'), 
                    os.path.join(py_eval_path, '1_full_sequence_dataset'), os.path.join(py_eval_path, '0_data_mat')]
        check =self.check_existence(paths)
        if check:
            self.clean(os.path.join(py_eval_path, '2_Cut_Data'), 2)
            self.clean(os.path.join(py_eval_path, '5_inferred'), 2)
            self.clean(os.path.join(py_eval_path, '6_results'), 1)
            self.clean(os.path.join(py_eval_path, '1_full_sequence_dataset'), 1)
            self.clean(os.path.join(py_eval_path, '0_data_mat'), 1)