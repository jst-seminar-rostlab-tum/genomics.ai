export default async function queryMyInstitutions() {
  return [
    {
      id: 1,
      name: 'Helmholtz Institute',
      country: 'Germany',
      profilePictureURL: 'https://www.hzdr.de/db/Pic?pOid=55058',
      backgroundPictureURL: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAAXNSR0IArs4c6QAAAA1JREFUGFdjsC9c9x8ABK8CXkzrf1cAAAAASUVORK5CYII=',
      adminIds: [],
    },
    {
      id: 2,
      name: 'Technische Universität München',
      country: 'Germany',
      profilePictureURL: 'https://scalings.eu/wp-content/uploads/2019/11/tum-logo.png',
      backgroundPictureURL: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAAXNSR0IArs4c6QAAAA1JREFUGFdjsC9c9x8ABK8CXkzrf1cAAAAASUVORK5CYII=',
      adminIds: [1],
    },
    {
      id: 3,
      name: 'Rostlab',
      country: 'Germany',
      profilePictureURL: 'https://avatars.githubusercontent.com/u/4093405?s=200&v=4',
      backgroundPictureURL: 'data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAAXNSR0IArs4c6QAAAA1JREFUGFdjsC9c9x8ABK8CXkzrf1cAAAAASUVORK5CYII=',
      adminIds: [],
    },
  ];
}
