# Interface for the models
class model_interface:

    def set_config(new_config):
        pass


    def get_from_config(key):
        pass


    # python3.9 scVI.py --input data/pancreas_normalized.h5ad -t -q


    def setup():
        """
        Set up the warnings filter and the figure parameters
        :return:
        """
        pass


    def pre_process_data():
        """
        loads the source_adata/target_adata and removes the sparsity
        :return source_adata, target_adata:
        """
        pass




    def create_model(source_adata, target_adata):
        """
        if there is already a pretrained model, nothing happens otherwise a new one will be trained
        :param source_adata:
        :param target_adata:
        :return model:
        """
        pass


    def setup_anndata(anndata):
        """
        Just because it's prettier that way
        :param anndata:
        :return:
        """
        pass

    def get_model(anndata):
        """
        Just because it's prettier that way
        :param anndata:
        :return model:
        """
        


    def compute_latent(model, adata):
        """
        computes the latent of a model with specific adata
        :param model:
        :param adata:
        :return latent:
        """
        pass


    def compute_query(anndata):
        """
        trains the model on a query and saves the result
        :param anndata:
        :return model:
        """


    def compute_full_latent(source_adata, target_adata, model):
        """
        basically just takes to datasets, concatenates them and then computes the latent and saves the result
        :param source_adata:
        :param target_adata:
        :param model:
        :return:
        """


    def compute(new_config):
        """
        computes the model
        :param config:
        :return model_attributes:
        """

