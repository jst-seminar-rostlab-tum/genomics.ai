import axiosInstance from './axiosInstance';

const MODEL = 'models';

const ModelService = {
  getModeles: async () => {
    const { data } = await axiosInstance.get(`/${MODEL}`);
    return data;
  },

  getModel: async (id) => {
    const { data } = await axiosInstance.get(`/model/${id}`);
    return data;
  },
};

export default ModelService;
