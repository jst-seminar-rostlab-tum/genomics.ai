import axiosInstance from './axiosInstance';
import { startOrContinueUpload } from './UploadLogic';
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

  startOrContinueProjectUpload: async (
    selectedFile,
    submissionProgress,
    setSubmissionProgress,
    projectData,
  ) => startOrContinueUpload(selectedFile,
    submissionProgress,
    setSubmissionProgress,
    projectData),

  // Temporary solution for search to use backend data while teampage uses mock
  getTeamProjects: async (teamId) => MockProjectService.getTeamProjects(teamId),
};

export default ProjectService;
