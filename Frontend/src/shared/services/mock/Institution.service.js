import helmholtz from 'assets/helmholtz-logo.jpg';
import tum from 'assets/tum-logo.png';
import MemberService from '../Member.service';
import ProfileService from '../Profile.service';
import TeamService from '../Team.service';

const mockLeftIds = [];
let runningId = 3;

let mockInstitutions = [
  {
    id: '1',
    name: 'Helmholtz Institute',
    country: 'Germany',
    description: 'The primary objective of Helmholtz Zentrum München in its role as German Center for Environmental Health is to improve the health of each individual, but also the health of society as a whole (public health). \n\nTo achieve this goal, the Center conducts research on the interaction of genetics, environmental factors and lifestyle and applies the insights gained from this research to develop individualized strategies for the diagnosis, treatment and prevention of serious common diseases such as diabetes mellitus, allergies and lung diseases.',
    backgroundPictureURL: 'https://upload.wikimedia.org/wikipedia/commons/thumb/2/22/Helmholtz_Zentrum_M%C3%BCnchen.jpg/2560px-Helmholtz_Zentrum_M%C3%BCnchen.jpg',
    adminIds: ['626bdb1ed76c8b968a50f833'],
    memberIds: [2],
    avatarUrl: helmholtz,
  },
  {
    id: '2',
    name: 'Technische Universität München',
    country: 'Germany',
    description: 'The Technical University of Munich (TUM) is one of Europe’s top universities. It is committed to excellence in research and teaching, interdisciplinary education and the active promotion of promising young scientists. The university also forges strong links with companies and scientific institutions across the world. \nTUM was one of the first universities in Germany to be named a University of Excellence. Moreover, TUM regularly ranks among the best European universities in international rankings.',
    backgroundPictureURL: 'https://www.in.tum.de/fileadmin/_processed_/5/5/csm_2006_1015Bild0136_9dc504e910.jpg',
    avatarUrl: tum,
    adminIds: ['626bdb1ed76c8b968a50f833'],
    memberIds: [1, 2, 3, 4, 5],
  },
  {
    id: '3',
    name: 'Rostlab',
    country: 'Germany',
    description: 'The lab\'s research is driven by a conviction that protein and DNA sequences encode a significant core of information about the ultimate structure and function of genetic material and its gene products. \nResearch goals of the lab involve using protein and DNA sequences along with evolutionary information to predict a protein\'s: overall function, interaction partners, secondary structure, disordered regions, subcellular localization, membrane spanning protein structure, intra-chain residue contacts, cell cycle control, and domain boundaries. \nAnother significant research focus is to improve the effectiveness and efficiency of structural genomics projects\' ability to determine the structures of proteins on a large scale.',
    avatarUrl: 'https://avatars.githubusercontent.com/u/4093405?s=200&v=4',
    adminIds: [1],
    memberIds: ['626bdb1ed76c8b968a50f833', 2, 3, 4, 5],
  },
];

const InstitutionService = {
  async createInstitution(name, description) {
    // fake effect
    await new Promise((resolve) => setTimeout(resolve, 1000));
    const user = await ProfileService.getProfile();
    runningId += 1;
    const newInstitution = {
      id: runningId.toString(),
      name,
      country: 'Germany',
      description,
      avatarUrl: null,
      backgroundPictureURL: null,
      adminIds: [user.id], // TODO: make sure that the backend puts my user ID here
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
    let preparedInstitution = mockInstitutions.map((institution) => ({
      ...institution,
      profilePictureURL: institution.avatarUrl,
    }));
    if (params.keyword) {
      preparedInstitution = preparedInstitution.filter(
        (team) => team.name.toLowerCase().includes(params.keyword.toLowerCase()),
      );
    }

    if (!params.sortBy || params.sortBy === 'name') {
      preparedInstitution = preparedInstitution.sort((a, b) => (`${a.name}`).localeCompare(b.name));
    }
    return preparedInstitution;
  },

  getTeamsOfInstitutionById: async (id) => TeamService.getInstitutionTeams(id)
  ,

};

export default InstitutionService;
