import MockInstitutionService from './mock/Institution.service';
import axiosInstance from './axiosInstance';
import { enhanceMember } from './Member.service';
import ProfileService from './Profile.service';

const MOCK_INSTUTITIONS = false;

const MODEL = 'institutions';

function enhanceInstitution(institution) {
  return { ...institution, id: institution._id };
}

const InstitutionService = MOCK_INSTUTITIONS ? MockInstitutionService : {
  async createInstitution(name, country) {
    const { data } = await axiosInstance.post('/institutions', {
      name,
      country,
    });
    return enhanceInstitution(data);
  },

  async getMyInstitutions() {
    const user = await ProfileService.getProfile();
    let { data } = await axiosInstance.get(`/users/${user.id}/institutions`);
    data = data.map(enhanceInstitution);
    return data;
  },

  async getMyAdminInstitutions() {
    const myInstitutions = await InstitutionService.getMyInstitutions();
    const user = await ProfileService.getProfile();
    return myInstitutions.filter((i) => i.adminIds.includes(user.id));
  },

  async getInstitution(institutionId) {
    const { data } = await axiosInstance.get(`/institutions/${institutionId}`);
    return enhanceInstitution(data);
  },

  async leaveInstitution(institutionId) {
    const user = await ProfileService.getProfile();
    try {
      await axiosInstance.delete(`/institutions/${institutionId}/join`, { data: { userId: user.id } });
    } catch (e) {
      throw Error(e.response.data);
    }
  },

  async getMembers(institutionId) {
    const { data } = await axiosInstance.get(`/institutions/${institutionId}/members`);
    return data.map(enhanceMember);
  },

  getInstitutions: async (params) => {
    const preparedParams = { ...params };
    if (!preparedParams.sortBy) {
      preparedParams.sortBy = 'name';
    }
    const { data } = await axiosInstance.get(`/${MODEL}`, { params: preparedParams });
    return data.map(enhanceInstitution);
  },

  getInstitutionById: async (id) => {
    // change later with getInstitution(...)
    const { data } = await axiosInstance.get(`/${MODEL}/${id}`);
    return data;
  },

  getTeamsOfInstitutionById: async (id) => {
    const { data } = await axiosInstance.get(`/${MODEL}/${id}/teams`);
    return data;
  },

  async updateDetails(institutionId, details) {
    await axiosInstance.put(`/institutions/${institutionId}`, { description: details });
  },

  async inviteMember(institutionId, invitedMail) {
    await axiosInstance.put(`/institutions/${institutionId}/invite`, { email: invitedMail });
  },

  async removeMemberFromInstitution(institutionId, memberId) {
    try {
    await axiosInstance.delete(`/institutions/${institutionId}/members/${memberId}`);
    } catch (e) {
      throw Error(e.response.data);
    }
  },

  async makeInstitutionAdmin(institutionId, memberId) {
    await axiosInstance.put(`/institutions/${institutionId}/admin`, { userId: memberId });
  },

  async removeInstitutionAdmin(institutionId, memberId) {
    await axiosInstance.delete(`/institutions/${institutionId}/admins/${memberId}`);
  }
};

export default InstitutionService;
