import axiosInstance from './axiosInstance';
import ProjectMock from './mock/projects';
import { startOrContinueUpload } from './UploadLogic';

const MODEL = 'ownprojects';

const ProjectService = {
  getProjects: async () => {
    const { data } = await axiosInstance.get(`/${MODEL}`);
    // const data = await ProjectMock.getProjects();
    return data;
  },

  getProject: async (id) => {
  // const { data } = await axiosInstance.get(`/project/${id}`);
  // ProjectMock.getProject(id)

    const { data } = await axiosInstance.get(`/${MODEL}`);
    const project = await data.find((p) => String(p._id) === String(id));
    return project;
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
