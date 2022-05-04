import axiosInstance from './axiosInstance';

const ATLASES = 'atlases';
const ATLAS = 'atlas';

const AtlasService = {
  getAtlases: async () => {
    const { data } = await axiosInstance.get(`/${ATLASES}`);
    return data;
  },
  getAtlasById: async (id) => {
    const { data } = await axiosInstance.get(`/${ATLAS}/${id}`);
    return data;
  },
};

export default AtlasService;
