import MockTeamService from './mock/Team.service';
import axiosInstance from './axiosInstance';
import ProfileService from './Profile.service';

const MOCK_TEAMS = false;
const MODEL = 'teams';

function enhanceTeam(team) {
  return { ...team, id: team._id, name: team.title };
}

const TeamService = MOCK_TEAMS ? MockTeamService : {
  async leaveTeam(teamId) {
    const user = await ProfileService.getProfile();
    try {
      await axiosInstance.delete(`/teams/${teamId}/join`, { data: { userId: user.id } });
    } catch (e) {
      throw Error(e.response.data);
    }
  },

  async createTeam(name, description, institutionId) {
    const { data } = await axiosInstance.post('/teams', {
      title: name,
      description,
      institutionId,
      visibility: 'PRIVATE',
    });
    return enhanceTeam(data);
  },

  async removeMemberFromTeam(teamId, memberId) {
    const team = await TeamService.getTeam(teamId);
    team.memberIds = team.memberIds.filter((mId) => mId !== memberId);
    await axiosInstance.post(`/teams/${teamId}`, { data: team });
  },

  async getTeam(teamId) {
    const { data } = await axiosInstance.get(`/teams/${teamId}`);
    return enhanceTeam(data);
  },

  async changeTeamDescription(teamId, description) {
    const { data } = await axiosInstance.put(`/teams/${teamId}`, { description });
    return data;
  },

  async changeTeamVisibility(teamId, visibility) {
    const { data } = await axiosInstance.put(`/teams/${teamId}`, { visibility });
    return data;
  },

  async getInstitutionTeams(institutionId) {
    const { data } = await axiosInstance.get(`/institutions/${institutionId}/teams`);
    return data.map(enhanceTeam);
  },

  async getMembers(teamId) {
    const { data } = await axiosInstance.get(`/teams/${teamId}/members`);
    return data.map(enhanceTeam);
  },

  async getMyTeams() {
    const user = await ProfileService.getProfile();
    let { data } = await axiosInstance.get(`/users/${user.id}/teams`);
    data = data.map(enhanceTeam);
    return data;
  },

  getTeams: async (params) => {
    const { data } = await axiosInstance.get(`/${MODEL}`, { params });
    return data.map(enhanceTeam);
  },

  addProject: async (teamId, projectId) => {
    const { data } = await axiosInstance.put(`/${MODEL}/${teamId}/add_project`, { projectId });
    return data;
  },
};

export default TeamService;
