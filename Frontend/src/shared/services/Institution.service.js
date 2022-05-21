import MockInstitutionService from './mock/Institution.service';
import axiosInstance from './axiosInstance';
import { enhanceMember } from './Member.service';
import ProfileService from './Profile.service';

const MOCK_INSTUTITIONS = true;

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
    let { data } = await axiosInstance.get(`/user/${user.id}/institutions`);
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

  async getMembers(institutionId) {
    return []; // TODO: enable once exists
    // eslint-disable-next-line no-unreachable
    const { data } = await axiosInstance.get(`/institutions/${institutionId}/members`);
    return data.map(enhanceMember);
  },

  getInstitutions: async (params) => {
    const { data } = await axiosInstance.get(`/${MODEL}`, { params });
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
};

export default InstitutionService;
