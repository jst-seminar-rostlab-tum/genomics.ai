import axiosInstance from './axiosInstance';

const MODEL = 'teams';

const TeamService = {
  getTeams: async (params) => {
    const { data } = await axiosInstance.get(`/${MODEL}`, { params });
    return data;
  },
};

export default TeamService;
