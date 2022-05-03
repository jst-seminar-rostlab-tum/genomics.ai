import axiosInstance from './axiosInstance';

const MODELS = 'models';

const ModelsService = {
  getModels: async () => {
    const { data } = await axiosInstance.get(`/${MODELS}`);
    return data;
  },
};

export default ModelsService;
