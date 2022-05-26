import Avatars from 'components/Avatars';
import React from 'react';

const MemberAvatars = ({ members }) => members.length !== 0 && (
<Avatars
  items={members.slice(0, 3).map(({ firstName, lastName, avatarUrl }) => ({ src: avatarUrl, alt: `${firstName} ${lastName}` }))}
  displayMoreIcon={members.length > 3}
/>
);

export default MemberAvatars;
