import axiosInstance from './axiosInstance';
import ProjectMock from './mock/projects';

const MODEL = 'atlases';

const AtlasService = {
  getAtlases: async () => {
    const { data } = await axiosInstance.get(`/${MODEL}`);
    // const data = await ProjectMock.getAtlases();
    return data;
  },

  getAtlas: async (id) => {
    const { data } = await axiosInstance.get(`/atlas/${id}`);
    return data;
  },
};

export default AtlasService;
