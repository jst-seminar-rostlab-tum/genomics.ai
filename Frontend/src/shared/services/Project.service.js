import axiosInstance from './axiosInstance';

const MODEL = 'projects';

const ProjectService = {
  getProjects: async (params) => {
    const { data } = await axiosInstance.get(`/${MODEL}`, { params });
    return data;
  },
};

export default ProjectService;
