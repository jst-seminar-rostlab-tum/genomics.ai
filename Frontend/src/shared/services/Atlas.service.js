import axiosInstance from './axiosInstance';

const MODEL = 'atlases';

const AtlasService = {
  getAtlases: async () => {
    const { data } = await axiosInstance.get(`/${MODEL}`);
    return data;
  },

  getAtlas: async (id) => {
    const { data } = await axiosInstance.get(`/atlas/${id}`);
    return data;
  },
};

export default AtlasService;
