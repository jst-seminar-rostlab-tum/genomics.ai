import axiosInstance from './axiosInstance';
import MockProjectService from './mock/Project.service';

const MODEL = 'projects';
const MOCK_PROJECTS = false;

const ProjectService = MOCK_PROJECTS ? MockProjectService : {
  getProjects: async (params) => {
    const { data } = await axiosInstance.get(`/${MODEL}`, { params });
    return data;
  },
  // Temporary solution for search to use backend data while teampage uses mock
  getTeamProjects: async (teamId, forPart) => MockProjectService.getTeamProjects(teamId, forPart),
};

export default ProjectService;
