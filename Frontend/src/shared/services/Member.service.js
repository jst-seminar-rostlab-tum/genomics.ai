import axiosInstance from './axiosInstance';

export function enhanceMember(member) {
  return { ...member, id: member._id };
}

const MemberService = {
};

export default MemberService;
