import helmholtz from 'assets/helmholtz-logo.jpg';
import tum from 'assets/tum-logo.png';
import MemberService from './Member.service';
import ProfileService from './Profile.service';
import TeamService from './Team.service';

const mockLeftIds = [];
let runningId = 3;

let mockInstitutions = [
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
    avatarUrl: tum,
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
    const newInstitution = {
      id: runningId.toString(),
      name,
      country: 'Germany',
      description,
      avatarUrl: null,
      backgroundPictureURL: null,
      adminIds: [1], // TODO: make sure that the backend puts my user ID here
      memberIds: [],
    };
    mockInstitutions = [...mockInstitutions, newInstitution];
    return newInstitution;
  },

  async leaveInstitution(institution) {
    mockLeftIds.push(institution.id);
  },

  async getMyInstitutions() {
    return mockInstitutions.filter((institution) => !mockLeftIds.includes(institution.id));
  },

  async getMyAdminInstitutions() {
    const myInstitutions = await InstitutionService.getMyInstitutions();
    const user = await ProfileService.getProfile();
    return myInstitutions.filter((i) => i.adminIds.includes(user.id));
  },

  async getInstitution(id) {
    return mockInstitutions.find((institution) => institution.id === id);
  },

  async getMembers(institutionId) {
    const institution = await InstitutionService.getInstitution(institutionId);
    return Promise.all([
      ...institution.adminIds,
      ...institution.memberIds,
    ].map(MemberService.getMember));
  },

  async removeMemberFromInstitution(institutionId, memberId) {
    const institution = mockInstitutions.find((i) => i.id === institutionId);
    if (!institution) {
      throw new Error('Institution not found');
    }
    institution.memberIds = institution.memberIds.filter((id) => id !== memberId);
  },

  getInstitutionById: async (id) => InstitutionService.getInstitution(`${id}`),

  getInstitutions: async (params) => {
    let preparedInstitution = mockInstitutions;
    if (params.keyword) {
      preparedInstitution = preparedInstitution.filter(
        (team) => team.name.toLowerCase().includes(params.keyword.toLowerCase()),
      );
    }
    return preparedInstitution.sort((a, b) => (`${a.name}`).localeCompare(b.name));
  },

  getTeamsOfInstitutionById: async (id) => TeamService.getInstitutionTeams(+id)
  ,

};

export default InstitutionService;
