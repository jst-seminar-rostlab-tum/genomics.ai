import helmholtz from 'assets/helmholtz-logo.jpg';
import getMember from './members';

const mockLeftIds = [];
let runningId = 3;

const testInstitutions = [
  {
    id: '1',
    name: 'Helmholtz Institute',
    country: 'Germany',
    description: 'Test',
    backgroundPictureURL: 'https://upload.wikimedia.org/wikipedia/commons/thumb/2/22/Helmholtz_Zentrum_M%C3%BCnchen.jpg/2560px-Helmholtz_Zentrum_M%C3%BCnchen.jpg',
    adminIds: ['626bdb1ed76c8b968a50f833'],
    memberIds: [2],
    avatarUrl: helmholtz,
  },
  {
    id: '2',
    name: 'Technische Universität München',
    country: 'Germany',
    backgroundPictureURL: 'https://www.in.tum.de/fileadmin/_processed_/5/5/csm_2006_1015Bild0136_9dc504e910.jpg',
    avatarUrl: 'https://scalings.eu/wp-content/uploads/2019/11/tum-logo.png',
    adminIds: [1],
    memberIds: ['626bdb1ed76c8b968a50f833', 2, 3, 4, 5],
  },
  {
    id: '3',
    name: 'Rostlab',
    country: 'Germany',
    avatarUrl: 'https://avatars.githubusercontent.com/u/4093405?s=200&v=4',
    adminIds: [1],
    memberIds: ['626bdb1ed76c8b968a50f833', 2, 3, 4, 5],
  },
];

const InstitutionService = {
  async createInstitution(name, description) {
    // fake effect
    await new Promise((resolve) => setTimeout(resolve, 1000));
    runningId += 1;
    return {
      id: runningId,
      name,
      country: null,
      description,
      avatarUrl: null,
      backgroundPictureURL: null,
      adminIds: [1], // TODO: make sure that the backend puts my user ID here
      memberIds: [],
    };
  },

  async leaveInstitution(institution) {
    mockLeftIds.push(institution.id);
  },

  async getMyInstitutions() {
    return testInstitutions.filter((institution) => !mockLeftIds.includes(institution.id));
  },

  queryIsAdminInstitutions(userId) {
    return testInstitutions.filter((institution) => institution.adminIds.includes(userId));
  },

  async getInstitution(id) {
    return testInstitutions.find((institution) => institution.id === id);
  },

  getInstitutionMembers: async (institutionId) => {
    const institution = await InstitutionService.getInstitution(institutionId);
    return Promise.all([...institution.adminIds, ...institution.memberIds].map(getMember));
  },

  async removeMemberFromInstitution(institutionId, memberId) {
    const institution = testInstitutions.find((i) => i.id === institutionId);
    if (!institution) {
      throw new Error('Institution not found');
    }
    institution.memberIds = institution.memberIds.filter((id) => id !== memberId);
  },
};

export default InstitutionService;
