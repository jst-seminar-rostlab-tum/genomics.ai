import axiosInstance from './axiosInstance';
import { startOrContinueUpload } from './UploadLogic';

const MODEL = 'ownprojects';

const ProjectService = {
  getProjects: async () => {
    const { data } = await axiosInstance.get(`/${MODEL}`);
    return data;
  },

  getProject: async (id) => {
    // const { data } = await axiosInstance.get(`/project/${id}`);
    const { data } = await axiosInstance.get(`/${MODEL}`);
    const project = await data.find((p) => String(p._id) === String(id));
    return { ...project, location: './testData/test_file1.csv' };
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
