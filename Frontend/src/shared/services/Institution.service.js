import axiosInstance from './axiosInstance';

const MODEL = 'institutions';

const InstitutionService = {
  getInstitutions: async (params) => {
    const { data } = await axiosInstance.get(`/${MODEL}`, { params });
    return data;
  },

  getInstitutionById: async (id) => {
    const { data } = await axiosInstance.get(`/${MODEL}/${id}`);
    return data;
  },

  getTeamsOfInstitutionById: async (id) => {
    const { data } = await axiosInstance.get(`/${MODEL}/${id}/teams`);
    return data;
  },
};

export default InstitutionService;
