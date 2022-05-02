import axiosInstance from './axiosInstance';
import ProjectMock from './mock/projects';

const MODEL = 'models';

const ModelService = {
  getModels: async () => {
    const { data } = await axiosInstance.get(`/${MODEL}`);
    // const data = await ProjectMock.getModels();
    return data;
  },

  getModel: async (id) => {
    const { data } = await axiosInstance.get(`/model/${id}`);
    return data;
  },
};

export default ModelService;
