import axiosInstance from './axiosInstance';
import { startOrContinueUpload } from './UploadLogic';

const MODEL = 'projects';

const ProjectService = {
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

  deleteItem: async (id) => {
    return ProjectService.deleteItem(id);
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
};

export default ProjectService;
