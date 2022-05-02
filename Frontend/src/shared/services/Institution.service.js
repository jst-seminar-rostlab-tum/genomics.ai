import axiosInstance from './axiosInstance';

const MODEL = 'institutions';

const InstitutionService = {
  getInstitutions: async (params) => {
    const { data } = await axiosInstance.get(`/${MODEL}`, { params });
    return data;
  },
};

export default InstitutionService;
