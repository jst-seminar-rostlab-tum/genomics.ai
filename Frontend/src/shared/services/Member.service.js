// import axiosInstance from './axiosInstance';
import MockMemberService from './mock/Member.service';

const MOCK_MEMBERS = false;

export function enhanceMember(member) {
  return { ...member, id: member._id };
}

const MemberService = MOCK_MEMBERS ? MockMemberService : {
};

export default MemberService;
