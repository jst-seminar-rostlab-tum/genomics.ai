import axiosInstance from './axiosInstance';
import { startOrContinueUpload } from './UploadLogic';

const MODEL = 'projects';

const ProjectService = {
  getProjects: async () => {
    const { data } = await axiosInstance.get(`/${MODEL}`);
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
};

export default ProjectService;
