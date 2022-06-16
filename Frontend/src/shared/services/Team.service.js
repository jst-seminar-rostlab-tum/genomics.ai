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
    try {
      await axiosInstance.delete(`/teams/${teamId}/join`);
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

  async joinTeam(teamId) {
    const user = await ProfileService.getProfile();
    const { data: updatedTeam } = await axiosInstance.put(`/${MODEL}/${teamId}/join`, { userId: user.id });
    return updatedTeam;
  },

  async inviteMemberByEmail(teamId, email) {
    try {
      await axiosInstance.put(`/teams/${teamId}/invite`, { email });
    } catch (e) {
      throw Error(e.response.data);
    }
  },

  async removeMemberFromTeam(teamId, memberId) {
    try {
      await axiosInstance.delete(`/teams/${teamId}/members/${memberId}`);
    } catch (e) {
      throw Error(e.response.data);
    }
  },

  async makeTeamAdmin(teamId, userId) {
    try {
      await axiosInstance.put(`/teams/${teamId}/admin`, { userId });
    } catch (e) {
      throw Error(e.response.data);
    }
  },

  async removeTeamAdmin(teamId, userId) {
    try {
      await axiosInstance.delete(`/teams/${teamId}/admins/${userId}`);
    } catch (e) {
      throw Error(e.response.data);
    }
  },

  async removeProjectFromTeam(teamId, projectId) {
    try {
      await axiosInstance.delete(`/teams/${teamId}/projects/${projectId}`);
    } catch (e) {
      throw Error(e.response.data);
    }
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
    const preparedParams = { ...params };
    if (!preparedParams.sortBy || preparedParams.sortBy === 'name') {
      preparedParams.sortBy = 'title';
    }
    const { data } = await axiosInstance.get(`/${MODEL}`, { params: preparedParams });

    return data.map(enhanceTeam);
  },

  addProject: async (teamId, projectId) => {
    const { data } = await axiosInstance.put(`/${MODEL}/${teamId}/add_project`, { projectId });
    return data;
  },
};

export default TeamService;
