import axiosInstance from './axiosInstance';
import MockProjectService from './mock/Project.service';

const MODEL = 'projects';
const MOCK_PROJECTS = false;

const ProjectService = MOCK_PROJECTS ? MockProjectService : {
  getProjects: async (params) => {
    const { data } = await axiosInstance.get(`/${MODEL}`, { params });
    return data;
  },

  getOwnProjects: async () => {
    const { data } = await axiosInstance.get('/ownprojects');
    return data;
  },

  getProject: async (id) => {
    const { data } = await axiosInstance.get(`/project/${id}`);
    return data;
  },

  createProject: async (projectName, atlasId, modelId, fileName) => {
    const { data } = await axiosInstance.post('/file_upload/start_upload', {
      projectName, atlasId, modelId, fileName,
    });
    console.log("The value that is passed is:" + JSON.stringify(data));
    return data;
  },

  deleteProject: async (id) => {
    await axiosInstance.delete(`/project/${id}`);
  },

  getDeletedProjects: async () => {
    const { data } = await axiosInstance.get('/deletedprojects');
    return data;
  },

  restoreProject: async (id) => {
    await axiosInstance.post(`/deletedprojects/${id}/restore`);
  },

  getTeamProjects: async (teamId) => {
    const { data } = await axiosInstance.get(`/teams/${teamId}/projects`);
    return data;
  }
};

export default ProjectService;
