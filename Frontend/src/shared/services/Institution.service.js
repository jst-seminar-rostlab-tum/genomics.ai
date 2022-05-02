import axiosInstance from './axiosInstance';
import { enhanceMember } from './Member.service';

function enhanceInstitution(institution) {
  return { ...institution, id: institution._id };
}

const InstitutionService = {
  getMyInstitutions: async () => {
    // const user = await ProfileService.getProfile();
    // TODO: change to /user/${user.id}/institutions
    let { data } = await axiosInstance.get('/institutions');
    data = data.map(enhanceInstitution);
    return data;
  },
  getInstitution: async (institutionId) => {
    const { data } = await axiosInstance.get(`/institutions/${institutionId}`);
    return enhanceInstitution(data);
  },
  getInstitutionMembers: async (institutionId) => {
    return []; // TODO: enable once exists
    // eslint-disable-next-line no-unreachable
    const { data } = await axiosInstance.get(`/institutions/${institutionId}/members`);
    return data.map(enhanceMember);
  },
};

export default InstitutionService;
