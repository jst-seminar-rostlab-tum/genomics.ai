import axiosInstance from './axiosInstance';

const MODEL = 'teams';

const TeamService = {
  addProject: async (teamId, projectId) => {
    const { data } = await axiosInstance.put(`/${MODEL}/${teamId}/add_project`, { projectId });
    return data;
  },
};

export default TeamService;
