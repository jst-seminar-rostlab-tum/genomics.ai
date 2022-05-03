import axiosInstance from './axiosInstance';
import ProjectMock from './mock/projects';

const MODEL = 'atlases';

const AtlasService = {
  getAtlases: async () => {
    const { data } = await axiosInstance.get(`/${MODEL}`);
    return data?.length ? data : ProjectMock.getAtlases();
  },

  getAtlas: async (id) => {
    // const { data } = await axiosInstance.get(`/atlas/${id}`);

    const { data } = await axiosInstance.get(`/${MODEL}`);
    const atlas = data.find((p) => String(p._id) === String(id));
    return atlas || ProjectMock.getAtlas(id);
  },
};

export default AtlasService;
