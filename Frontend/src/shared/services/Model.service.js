import axiosInstance from './axiosInstance';
import ProjectMock from './mock/projects';

const MODEL = 'models';

const ModelService = {
  getModels: async () => {
    const { data } = await axiosInstance.get(`/${MODEL}`);
    return data?.length ? data : ProjectMock.getModels();
  },

  getModel: async (id) => {
    // const { data } = await axiosInstance.get(`/model/${id}`);

    const { data } = await axiosInstance.get(`/${MODEL}`);
    const model = data.find((p) => String(p._id) === String(id));
    return model || ProjectMock.getModel(id);
  },
};

export default ModelService;
