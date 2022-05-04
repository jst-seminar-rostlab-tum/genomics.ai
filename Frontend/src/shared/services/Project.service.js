import axiosInstance from './axiosInstance';
import ProjectMock from './mock/projects';
import { startOrContinueUpload } from './UploadLogic';

const MODEL = 'ownprojects';

const ProjectService = {
  getProjects: async () => {
    const { data } = await axiosInstance.get(`/${MODEL}`);
    const mockData = await ProjectMock.getProjects();
    return [...data, ...mockData];
  },

  getProject: async (id) => {
  // const { data } = await axiosInstance.get(`/project/${id}`);

    const { data } = await axiosInstance.get(`/${MODEL}`);
    const project = data.find((p) => String(p._id) === String(id));
    return project || ProjectMock.getProject(id);
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
