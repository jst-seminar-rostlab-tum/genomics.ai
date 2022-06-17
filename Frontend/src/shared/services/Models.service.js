import axiosInstance from './axiosInstance';

const MODELS = 'models';
const MODEL = 'model';

const ModelsService = {
  getModels: async () => {
    const { data } = await axiosInstance.get(`/${MODELS}`);
    return data;
  },
  getModelById: async (id) => {
    const { data } = await axiosInstance.get(`${MODEL}/${id}`);
    return data;
  }
};

export default ModelsService;
