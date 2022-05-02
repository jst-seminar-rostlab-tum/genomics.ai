import MockInstitutionService from './mock/Institution.service';
import axiosInstance from './axiosInstance';
import { enhanceMember } from './Member.service';
import ProfileService from './Profile.service';

const MOCK_INSTUTITIONS = true;

function enhanceInstitution(institution) {
  return { ...institution, id: institution._id };
}

const InstitutionService = MOCK_INSTUTITIONS ? MockInstitutionService : {
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
};

export default InstitutionService;
