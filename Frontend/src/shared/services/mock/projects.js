export default async function getProject(id) {
  await new Promise((resolve) => setTimeout(resolve, 1000));
  return {
    async json() {
      return {
        id,
        name: 'Demo Project',
        resultURL: './testData/test_file1.csv',
      };
    },
  };
}
